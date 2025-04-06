using LinearAlgebra
using PrettyTables
using Distributions
using Clustering
using DataFrames
using Statistics
using Random
using Plots
using CSV
using SparseArrays  # Agregado para matrices dispersas

# ==============================================================================
# Lectura de datos (sin cambios)
# ==============================================================================
lines = DataFrame(CSV.File(joinpath(@__DIR__, "lines.csv")))
compensation = DataFrame(CSV.File(joinpath(@__DIR__, "compensation.csv")))
trafos = DataFrame(CSV.File(joinpath(@__DIR__, "trafos.csv")))
loads = DataFrame(CSV.File(joinpath(@__DIR__, "loads.csv")))
gen = DataFrame(CSV.File(joinpath(@__DIR__, "gen.csv")))

# ==============================================================================
# Función optimizada para calcular Ynn y Yns
# ==============================================================================
function calcular_Ynn_Yns(lines, comp, trafos)
    """
    Calcula Ynn y Yns para cada nodo y subestación.
    Entradas : 
              Lines : Dataframe con la información de las líneas.
              Trafos : Dataframe con la información de las subestaciones.
    Salidas : 
              Ynn : Matriz de admisibilidad de nodos.
              Yns : Matriz de admisibilidad de subestaciones.
              Ybus : Matriz de admisibilidad de nodos y subestaciones.
    """
    num_nodes = max(maximum(lines.from), maximum(lines.to))
    # Usar matriz dispersa en lugar de densa
    Ybus = spzeros(ComplexF64, num_nodes, num_nodes)
    
    s = 788
    sbase = 10^5

    # Precalcular índices para líneas
    line_idx = collect(1:nrow(lines))
    for k in line_idx
        n1, n2 = lines.from[k], lines.to[k]
        yL = 1/(lines.r_pu[k] + lines.x_pu[k]*1im)
        Ybus[n1,n1] += yL
        Ybus[n1,n2] -= yL
        Ybus[n2,n1] -= yL
        Ybus[n2,n2] += yL
    end

    # Vectorizar compensaciones
    for k in 1:nrow(comp)
        Ybus[comp.node[k], comp.node[k]] += comp.q_kVAr[k]/sbase*1im
    end

    # Optimización de transformadores
    for k in 1:nrow(trafos)
        n1, n2 = trafos.from[k], trafos.to[k]
        t = 1 + (trafos.position[k]*((trafos.t_max[k]-trafos.t_min[k])/(trafos.steps[k]-1)))
        yL = 1/(trafos.r_pu[k] + trafos.x_pu[k]*1im)
        Ybus[n1,n1] += (t^2 - t)*yL + t*yL
        Ybus[n1,n2] -= t*yL
        Ybus[n2,n1] -= t*yL
        Ybus[n2,n2] += (1-t)*yL + t*yL
    end

    # Separación eficiente
    idx = setdiff(1:num_nodes, s)
    Ynn = Ybus[idx, idx]
    Yns = Ybus[idx, s]
    return Ynn, Yns, Ybus
end

# ==============================================================================
#  Flujo estático optimizado
# ==============================================================================
function Flujo_punto_fijo(lines, compensation, trafos, loads, gen)
    """
    Función que calcula el flujo estático optimizado de la red.
    Entradas : 
            lines : Dataframe con la información de las líneas.
            compensation : Dataframe con la información de la compensación.
            trafos : Dataframe con la información de los transformadores.
            loads : Dataframe con la información de las cargas.
            gen : Dataframe con la información de las generaciones.
    Salidas :   
            Vns : Vector con los valores de tensión en los nodos.
    """
    num_nodes = max(maximum(lines.from), maximum(lines.to))
    sbase = 10^5
    
    Ynn, Yns, Ybus = calcular_Ynn_Yns(lines, compensation, trafos)
    iYnn = factorize(Ynn)  # Factorización LU de la matriz dispersa
    
    Vn = ones(ComplexF64, num_nodes-1)
    Vs = 1.0 + 0.0im
    
    # Precalcular Sn
    Sn = zeros(ComplexF64, num_nodes-1)
    # Se agrega la generación al vector de potencias inyectadas
    for k in 1:nrow(gen)
        n1 = gen.node[k]
        Sn[n1-1] += gen.Pg[k]/10^6 + (gen.Qg[k]/10^6)*1im
    end
    # Se agrega la carga al vector de potencias inyectadas
    for k in 1:nrow(loads)
        n1 = loads.nodo[k]
        Sn[n1-1] += -loads.p_kW[k]/sbase - (loads.q_kVAr[k]/sbase)*1im
    end
    
    # Iteración con criterio de convergencia
    max_iter = 10
    tol = 1e-6
    for i in 1:max_iter
        Vn_ant = copy(Vn)
        # Calcular el lado derecho como vector denso
        rhs = conj.(Sn ./ Vn) .- Yns * Vs
        # Resolver el sistema usando la factorización (convierte a denso si es necesario)
        Vn = iYnn \ Vector(rhs)  # Convertimos explícitamente a vector denso
        if norm(Vn - Vn_ant) < tol
            break
        end
    end
    Vns = zeros(num_nodes)*1im
    s = 788
    # Se agrega el nodo slack
    Vns[1: s-1] = Vn[1: s-1]
    Vns[s] = Vs
    Vns[s+1: end] = Vn[s:end]
    
    df = DataFrame(Nodos = 1:num_nodes, Vn = Vns)
    pretty_table(df)
    return Vns
end

# ==============================================================================
# Creando los perfiles de demanda a partir de un random walk
# ==============================================================================
function random_walk_complex(minutes, nodes, distribution, params)
    """
    Simula un random walk para generar perfiles de demanda para cada nodo.
    Entradas : 
            minutes : Valor int que representa el tiempo de simulación en minutos.
            nodes : Valor int que representa el número de nodos en la red.
            distribution : Valor string que representa la distribución de probabilidad de cada paso del random walk.
            params : Valor array que contiene los parámetros de la distribución.
    Salidas :
            matriz_rw : Matriz de tamaño (minutes, nodes) que contiene la información de la demanda de cada nodo 
            en los minutos establecidos
    """
    valor_inicial = [(rand() + rand()*1im) * 0.35 for _ in 1:nodes]  # Un vector de valores iniciales entre [0, 0.2] + i[0, 0.2]
    
    # Inicializar la matriz con el valor inicial en todas las posiciones
    matriz_rw = Array{ComplexF64}(undef, minutes, nodes)
    matriz_rw[1, :] = valor_inicial  # Asignamos el vector de valores iniciales a la primera fila

    # Generar los pasos aleatorios para la parte real e imaginaria
    if distribution == "normal"
        mu = get(params, "mu", 0.0)
        sigma = get(params, "sigma", 1.0)
        dist = Normal(mu, sigma)
        pasos_real = rand(dist, minutes - 1, nodes)
        pasos_imag = rand(dist, minutes - 1, nodes) .* 1im
    elseif distribution == "uniforme"
        low = get(params, "low", -1.0)
        high = get(params, "high", 1.0)
        dist = Uniform(low, high)
        pasos_real = rand(dist, minutes - 1, nodes)
        pasos_imag = rand(dist, minutes - 1, nodes) .* 1im
    else
        error("Distribución no soportada. Usa 'normal' o 'uniforme'.")
    end

    # Combinar las partes real e imaginaria
    pasos_complex = pasos_real .+ pasos_imag

    # Llenar la matriz con el random walk
    for t in 2:minutes
        matriz_rw[t, :] = matriz_rw[t-1, :] + pasos_complex[t-1, :]
    end

    # Limitar la magnitud
    matriz_rw = map(x -> (abs(x) < 0.1 ? 0.1/abs(x)*x : (abs(x) > 4 ? 4/abs(x)*x : x)), matriz_rw)

    return matriz_rw
end

# ==============================================================================
# Flujo cuasidinámico optimizado
# ==============================================================================
function flujo_cuasidinamico(Ynn, Yns, rw, gen)
    """
    Función que simula el flujo cuasidinamico en una red dada la generación y demanda a lo largo de un periodo de tiempo.
    Entradas : 
              Ynn : Matriz de admisibilidad de nodos.
              Yns : Matriz de admisibilidad de subestaciones y nodos.
              rw : Matriz de random walk de tamaño (minutes, nodes) que contiene la información dada por el random walk.
              gen : Matriz de generación de tamaño (minutes, nodes) que contiene la información dada por la generación.
    Salidas :  
              Vns_df : Dataframe que contiene la información de flujo de carga en cada nodo a lo largo del tiempo.
    """
    num_nodes = max(maximum(lines.from), maximum(lines.to))
    sbase = 10^5
    iYnn = factorize(Ynn)
    
    Vn_s = Matrix{ComplexF64}(undef, num_nodes, nrow(rw))
    Sn = zeros(ComplexF64, num_nodes-1)
    gen_dict = Dict(gen.node[i] => i for i in 1:nrow(gen))
    
    @views for j in 1:nrow(rw)
        Vn = ones(ComplexF64, num_nodes-1)
        Vs = 1.0 + 0.0im
        
        for k in 1:ncol(rw)
            nodo = parse(Int, names(rw)[k][2:end]) - 1
            if haskey(gen_dict, nodo)
                n = gen_dict[nodo]
                Sn[nodo] = (gen.Pg[n] - real(rw[j,k]))/sbase + ((gen.Qg[n] - imag(rw[j,k]))/sbase)*1im
            else
                Sn[nodo] = -rw[j,k]/sbase
            end
        end
        
        for r in 1:4
            # Convertir el lado derecho a vector denso explícitamente
            rhs = Vector(conj.(Sn ./ Vn) .- Yns * Vs)
            Vn = iYnn \ rhs
        end
        Vs = 1 + 0*1im
        Vns = zeros(num_nodes)*1im
        s = 788
        # Se agrega el nodo slack
        Vns[1: s-1] = Vn[1: s-1]
        Vns[s] = Vs
        Vns[s+1: end] = Vn[s:end]
        Vn_s[:,j] = Vns
    end 
    Vns_df = DataFrame(Vn_s, :auto)
    return Vns_df
end

# ==============================================================================
# Ajuste de taps optimizado
# ==============================================================================
function Ajuste_aut_taps(lines, comp, trafos, rw, gen)
    """
    Función para el ajuste de taps optimizado para cada transformador(taps en el LV) de la red.
    Entradas : 
            lines : Dataframe que contiene la información de los transformadores de la red.
            comp : Dataframe que contiene la información de los nodos de la red.
            trafos : Dataframe que contiene la información de los transformadores de la red.
            rw : Dataframe que contiene la información de los perfiles de demanda generados por el random walk.
            gen : Dataframe que contiene la información de los generadores de la red.
    Salidas : 
            taps_df : Dataframe que contiene la información de los taps optimizados para cada transformador de la red.
            Vns_5min : Dataframe que contiene la información de los valores de tensión en cada nodo de 
            la red para cada intervalo de 5 minutos.
    """
    Ynn, Yns, _ = calcular_Ynn_Yns(lines, comp, trafos)
    Vns_base = flujo_cuasidinamico(Ynn, Yns, rw, gen)
    n = max(maximum(lines.to), maximum(lines.from))-1
    Vns_5min = DataFrame()
    taps_df = DataFrame()
    trafo_dict = Dict(trafos.to[i] => i for i in 1:nrow(trafos))
    int = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    for i in int
        trafos_c = deepcopy(trafos)
        println("Iter: ", i)
        for j in trafos_c.to
            display(j)
            while true
                ind = trafo_dict[j]
                V = abs(Vns_base[j,i])
                steps = trafos.steps[ind]
                pos = trafos_c.position[ind]
                println("Vns: " , V)
                println("Pos_tap: ", pos)
                if V > 1.05 && pos > -(steps-1)/2
                    trafos_c.position[ind] -= 1
                elseif V < 0.95 && pos < (steps-1)/2
                    trafos_c.position[ind] += 1
                else
                    break
                end
                    Ynn, Yns, _ = calcular_Ynn_Yns(lines, comp, trafos_c)
                    Vns_base = flujo_cuasidinamico(Ynn, Yns, rw, gen)
            end
            ind = trafo_dict[j]
            V = abs(Vns_base[j,i])
            pos = trafos_c.position[ind]
           
        end
        taps_df[!, "taps_int_min_$i"] = deepcopy(trafos_c.position)
        Vns_5min[!, "Vns_min$i"] = deepcopy(Vns_base[:,i])
    end
    return taps_df, Vns_5min
end

# ==============================================================================
# Ejecución del punto fijo
# ==============================================================================
Vns = Flujo_punto_fijo(lines, compensation, trafos, loads, gen)


# ==============================================================================
# Creación del Random Walk
# ==============================================================================
Ynn, Yns, Ybus = calcular_Ynn_Yns(lines, compensation, trafos)
minutos, nodes = 60, nrow(loads)
rw = random_walk_complex(minutos, nodes, "normal", Dict("mu" => 0.05, "sigma" => 0.005))
rw_df = DataFrame(rw, ["N$(loads.nodo[i])" for i in 1:nrow(loads)])
CSV.write(joinpath(@__DIR__, "random_walk.csv"), rw_df)

# ==============================================================================
# Calculo de pérdidas 
# ==============================================================================
sbase = 10^5
Perd = Vns'*Ybus*Vns
display(Perd*sbase)

# ==============================================================================
# Ejecución del flujo cuasidinamico caso base
# ==============================================================================
Vns_base = flujo_cuasidinamico(Ynn, Yns, rw_df, gen)
CSV.write(joinpath(@__DIR__, "V_fcd_base.csv"), abs.(Vns_base))

# ==============================================================================
# Se ejecuta el ajuste automatico de taps
# ==============================================================================
@time taps_df, Vns_base = Ajuste_aut_taps(lines, compensation, trafos, rw_df, gen)
CSV.write(joinpath(@__DIR__, "V_fcd_taps.csv"), abs.(Vns_base))
CSV.write(joinpath(@__DIR__, "taps_corregidos.csv"), taps_df)
tabla1 = pretty_table(abs.(Vns_base))
display(tabla1)
tabla2 = pretty_table(taps_df)
display(tabla2)