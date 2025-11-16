module GetExclusiveKmers

using FASTX

export GetExclusiveKmers

function loadStringSequences(
    file::String
)::Vector{String}

    sequences = String[]
    for record in open(FASTAReader, file)
        push!(sequences, sequence(String, record))
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
    base = UInt64(257)

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




function getOccursin_rolling_hash(
    sequence::String,
    kmer_hash_map::Dict{UInt64,String},
    k_len::Int)::Dict{UInt64,String}

    seq_len = length(sequence)

    update_hashmap = Dict{UInt64,String}()

    if seq_len < k_len
        error("K-mer size must be greater than the inputed sequence")
    end

    base = UInt64(257)

    # Calculate base^(k_len-1) for rolling hash
    power = UInt64(1)
    for i in 1:(k_len-1)
        power *= base
    end

    # Calculate initial hash
    current_hash = UInt64(0)
    @inbounds for i in 1:k_len
        current_hash = current_hash * base + UInt64(sequence[i])
    end

    # Check first k-mer
    if haskey(kmer_hash_map, current_hash)
        update_hashmap[current_hash] = kmer_hash_map[current_hash]
    end

    # Roll through sequence
    @inbounds for i in (k_len+1):seq_len
        # Rolling hash: remove leftmost char, add rightmost char
        current_hash = current_hash - UInt64(sequence[i-k_len]) * power
        current_hash = current_hash * base + UInt64(sequence[i])

        if haskey(kmer_hash_map, current_hash)
            update_hashmap[current_hash] = kmer_hash_map[current_hash]
        end
    end

    return update_hashmap

end



function julia_main()::Cint

    # ACTG


    # CTGA
    # TGAC
    # GACT


    seqteste::Vector{String} = ["ACTGACT", "TCTGACT"]
    k::Int = 9
    numero_esperado = 4^k

    sketch::Dict{String,Int} = gerar_kmers(k)

    println("Total de k-mers esperados: $numero_esperado")
    println("Total de k-mers gerados: ", length(sketch))

    # Pressuposto - A população de sequencias possuem todos k-mers contidos
    # A cada iteração eu tenho que filtrar meu conjunto
    # Até acabar com as sequencias assim o que sobrar será o conjunto de k-mers exclusivo

    kmer_hash_map = Dict{UInt64,String}()

    for kmer in keys(sketch)
        h = compute_hash(kmer)
        if !haskey(kmer_hash_map, h)
            kmer_hash_map[h] = kmer
        end
    end

    # @show kmer_hash_map

    variantDirPath = "/home/salipe/Desktop/datasets/mkpx/data/"

    variantDirs::Vector{String} = readdir(variantDirPath)

    @inbounds for v in eachindex(variantDirs)
        variant::String = variantDirs[v]
        println("Processing $variant")
        # var_hash = kmer_hash_map


        sequences::Vector{String} = loadStringSequences("$variantDirPath/$variant/$variant.fasta")

        for seq in sequences
            kmer_hash_map = getOccursin_rolling_hash(seq, kmer_hash_map, k)
        end

        @show length(kmer_hash_map)
    end


    return 0

end

end


GetExclusiveKmers.julia_main()