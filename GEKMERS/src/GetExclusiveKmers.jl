module GetExclusiveKmers

using FASTX, ArgParse

export GetExclusiveKmers

function loadStringSequences(
    file::String
)::Vector{String}

    sequences = String[]
    for record in open(FASTAReader, file)
        push!(sequences, uppercase(sequence(String, record)))
    end
    return sequences
end


function loadCodeUnitsSequences(
    file::String
)::Vector{Base.CodeUnits}

    sequences = Vector{Base.CodeUnits}()
    for record in open(FASTAReader, file)
        push!(sequences, codeunits(sequence(String, record)))
    end
    return sequences
end

function compute_hash(s::String)::UInt64
    h = UInt64(0)
    base = UInt64(5)

    for char in s
        h = h * base + UInt64(char)
    end

    return h
end

function gerar_kmers(k::Int)

    nucleotides = ['A', 'C', 'T', 'G']

    # Inicializar o sketch como um dicionário vazio
    sketch = Dict{String,Int}()

    # Função recursiva para gerar todas as combinações
    function gerar_combinacoes(prefixo::String, tamanho_restante::Int)
        if tamanho_restante == 0
            # Adicionar ao sketch (inicializando com valor 0)
            sketch[prefixo] = 0
            return
        end

        # Para cada nucleotídeo, adicionar ao prefixo e continuar recursivamente
        for nucleotide in nucleotides
            gerar_combinacoes(prefixo * nucleotide, tamanho_restante - 1)
        end
    end

    gerar_combinacoes("", k)

    return sketch
end


function max_entropy(kmers::Dict{String,Int32})
    # Extrair valores e ordenar
    data = collect(values(kmers))
    total = sum(data)
    sort!(data, rev=true)

    # Normalizar dados
    normalized_data = data ./ total


    function entropy(probs::Vector{Float64})::Float64
        -sum(p * log(p) for p in probs if p > 0.0)
    end

    # Calcular curva de entropia 
    entropy_curve = map(1:length(normalized_data)-1) do s
        # Região A: probs[1:s], probabilidade total P_A
        p_a = sum(normalized_data[1:s])
        h_a = if p_a > 0.0
            p_a_data = normalized_data[1:s] ./ p_a
            entropy(p_a_data)
        else
            0.0
        end

        # Região B: probs[s+1:end], probabilidade total P_B
        p_b = sum(normalized_data[s+1:end])
        h_b = if p_b > 0.0
            p_b_data = normalized_data[s+1:end] ./ p_b
            entropy(p_b_data)
        else
            0.0
        end

        h_a + h_b
    end

    # Encontrar índice de máxima entropia
    max_entropy_idx = argmax(entropy_curve)
    threshold = max_entropy_idx
    frequency = data[threshold]

    return frequency
end


function rolling_hash_kmers(
    sequence::String,
    kmer_hash_map::Dict{UInt64,Tuple{String,Int32}},
    k_len::Int)::Dict{UInt64,Tuple{String,Int32}}

    seq_len = length(sequence)

    update_hashmap = Dict{UInt64,Tuple{String,Int32}}()

    if seq_len < k_len
        error("K-mer size must be greater than the inputed sequence")
    end
    if length(keys(kmer_hash_map)) == 0
        error("0 k-mers found!")
    end

    base = UInt64(5)

    power = UInt64(1)
    for i in 1:(k_len-1)
        power *= base
    end


    # Have to init the hash
    current_hash = UInt64(0)
    @inbounds for i in 1:k_len
        current_hash = current_hash * base + UInt64(sequence[i])
    end

    if haskey(kmer_hash_map, current_hash)
        update_hashmap[current_hash] = kmer_hash_map[current_hash]
        kmer, value = update_hashmap[current_hash]
        update_hashmap[current_hash] = (kmer, value + 1)
    end

    @inbounds for i in (k_len+1):seq_len
        current_hash = current_hash - UInt64(sequence[i-k_len]) * power
        current_hash = current_hash * base + UInt64(sequence[i])


        if haskey(kmer_hash_map, current_hash)
            if haskey(update_hashmap, current_hash)
                kmer, value = update_hashmap[current_hash]
                update_hashmap[current_hash] = (kmer, value + 1)
            else
                update_hashmap[current_hash] = kmer_hash_map[current_hash]
                kmer, value = update_hashmap[current_hash]
                update_hashmap[current_hash] = (kmer, value + 1)
            end
        end
    end

    return update_hashmap

end

function generate_kmers(k::Int)
    nucleotides = ['A', 'C', 'T', 'G']
    sketch = Dict{String,Int}()
    function gerar_combinacoes(prefixo::String, tamanho_restante::Int)
        if tamanho_restante == 0
            sketch[prefixo] = 0
            return
        end

        for nucleotide in nucleotides
            gerar_combinacoes(prefixo * nucleotide, tamanho_restante - 1)
        end
    end

    gerar_combinacoes("", k)

    return sketch
end



function get_exclusive_kmers(
    k_len::Int,
    dataset_path::String,
    use_gramep::Bool,
    use_entropy::Bool,
    reference::Union{Nothing,String})::Set{String}

    sketch::Dict{String,Int} = generate_kmers(k_len)

    kmer_hash_map = Dict{UInt64,Tuple{String,Int32}}()

    for kmer in keys(sketch)
        h = compute_hash(kmer)
        if !haskey(kmer_hash_map, h)
            kmer_hash_map[h] = (kmer, Int32(0))
        end
    end

    all_exclusive_kmers = Set{String}()

    variantDirs::Vector{String} = readdir(dataset_path)
    variant_kmers = Dict{String,Set{String}}()

    # Calcular os kmers da referencia
    if use_gramep
        reference_seq::String = loadStringSequences(reference)[1]
        reference_var_hash = rolling_hash_kmers(reference_seq, kmer_hash_map, k_len)

        reference_kmer_dict = Dict{String,Int32}()
        for kmer_freq in values(reference_var_hash)
            reference_kmer_dict[kmer_freq[1]] = kmer_freq[2]
        end

        freq = max_entropy(reference_kmer_dict)
        filter!(e -> e[2] >= freq, reference_kmer_dict)
    end


    @inbounds for v in eachindex(variantDirs)
        variant::String = variantDirs[v]
        println("Getting $variant k-mers")
        var_hash = kmer_hash_map

        sequences::Vector{String} = loadStringSequences("$dataset_path/$variant/$variant.fasta")
        for seq in sequences
            var_hash = rolling_hash_kmers(seq, var_hash, k_len)
        end

        # A busca das mais informativas tem que ser aqui em comparação com a referencia
        kmer_dict = Dict{String,Int32}()
        for kmer_freq in values(var_hash)
            kmer_dict[kmer_freq[1]] = kmer_freq[2]
        end

        if length(keys(kmer_dict)) <= 1
            error("Insufficient k-mers found!")
        else
            @info "Found $(length(keys(kmer_dict))) kmers for $variant"
        end

        freq = max_entropy(kmer_dict)
        filter!(e -> e[2] >= freq, kmer_dict)

        if use_gramep
            filter!(e -> !haskey(reference_kmer_dict, e[2]), kmer_dict)
        end

        variant_kmers[variant] = Set(keys(kmer_dict))

        @info "Selected $(length(keys(kmer_dict))) kmers for $variant"

        for kmer_values in values(var_hash)
            union!(all_exclusive_kmers, Set([kmer_values[1]]))
        end

    end

    # Measure all te intersections
    kmers_in_multiple_variants = Set{String}()
    variant_list = collect(keys(variant_kmers))

    @inbounds for i in eachindex(variant_list)
        for j in (i+1):length(variant_list)
            v1 = variant_list[i]
            v2 = variant_list[j]
            shared = intersect(variant_kmers[v1], variant_kmers[v2])
            union!(kmers_in_multiple_variants, shared)
        end
    end

    union_exclusive = setdiff(all_exclusive_kmers, kmers_in_multiple_variants)
    return union_exclusive
end



function julia_main()::Cint


    settings = ArgParseSettings(
        description="GEK-MERS - Get Exclusive K-mers",
        version="0.1",
        add_version=true
    )

    @add_arg_table! settings begin
        "-g", "--use-gramep"
        help = "Use gramep selection"
        action = :store_true
        "-e", "--use-entropy"
        help = "Use Max Entropy to select k-mers"
        action = :store_true
        "-k", "--k-len"
        help = "K-mer K value"
        required = true
        arg_type = Int
        "-d", "--directory"
        help = "Dataset directory path"
        required = true
        arg_type = String
        "-r", "--reference"
        help = "Dataset directory path"
        required = false
    end

    parsed_args = parse_args(settings)

    directory::String = parsed_args["directory"]
    k_len::Int = parsed_args["k-len"]
    use_gramep::Bool = parsed_args["use-gramep"]
    use_entropy::Bool = parsed_args["use-entropy"]
    reference::Union{Nothing,String} = parsed_args["reference"]  # optional, may be `nothing`

    numero_esperado = 4^k_len
    @info "Quantidade de K-mer máximo:" numero_esperado

    exclusive_kmers::Set{String} = get_exclusive_kmers(k_len, directory, use_gramep, use_entropy, reference)

    @info "Found $(length(exclusive_kmers)) kmers for the dataset"


    return 0

end

end


GetExclusiveKmers.julia_main()