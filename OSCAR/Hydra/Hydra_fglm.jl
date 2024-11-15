using Oscar
include("Hydra.jl")
include("Hydra_polynomial_model.jl")

function is_in_shape_position(gb_lex, P)
    in_shape_position = true
    with_ordering(P, lex(P)) do
        for poly in gb_lex
            if !is_univariate(poly)
                if !is_univariate(poly - leading_term(poly))
                    in_shape_position = false
                    break
                end
            end
        end
    end
    return in_shape_position
end

function extract_Hydra_gb(hydra)
    m = 2
    polys = generate_Hydra_polynomials_m_samples(hydra=hydra, m=m, info_level=0);
    polys = transform_Hydra_polynomial_system(hydra, polys, m);
    _, _, polys_downsized_subs = non_linear_variable_substitution_Hydra_polynomial_system(hydra,
                                                                                          polys, 
                                                                                          m; 
                                                                                          transformed=true,
                                                                                          info_level=0);
    P = parent(polys_downsized_subs[1])
    variables_subs = map(i -> "x_subs_i" * string(i), 1:2 * hydra.rounds_head - 2)
    P_subs, variables_subs = polynomial_ring(hydra.field, variables_subs, internal_ordering=:degrevlex)
    zero_vec = zero_matrix(P_subs, length(gens(P)) - length(variables_subs), 1)
    image = [vec(zero_vec[:, 1]); variables_subs]
    phi = hom(P, P_subs, image)
    polys_downsized_subs = map(phi, polys_downsized_subs)
    gb_subs = [polys_downsized_subs[2 + 1:hydra.rounds_head];
               polys_downsized_subs[2 + hydra.rounds_head + 1:2 * hydra.rounds_head + 2]
               ]
    return gb_subs
end

function bench_fglm(gb_DRL, P; print_time=true)
    t_1 = time()
    gb_lex = fglm(ideal(gb_DRL), destination_ordering=lex(P))
    t_2 = time()
    t_lex = t_2 - t_1
    if print_time
        println("Time needed for FGLM: ", t_lex, "s")
        println("Shape position: ", is_in_shape_position(gb_lex, P))
    end
end

function random_quadratic_polynomial(P)
    K = base_ring(P)
    variables = gens(P)
    quadratic_monoms = vec(map(prod -> prod[1] * prod[2], Iterators.product(variables, variables)))
    poly = [map(mon -> rand(K) * mon, quadratic_monoms);
            map(mon -> rand(K) * mon, variables);
            [rand(K)]
            ]
    return sum(poly)
end

function random_quadratic_DRL_gb(P; nr_thrds=1, print_time=true)
    n = length(gens(P))
    polys = map(i -> random_quadratic_polynomial(P), 1:n)
    t_1 = time()
    gb_drl = groebner_basis_f4(ideal(polys), nr_thrds=nr_thrds, info_level=0)
    t_2 = time()
    t_drl = t_2 - t_1
    if print_time
        println("Time needed for random DRL Gröbner basis: ", t_drl, "s")
    end
    return gb_drl
end


nr_thrds = 48
K = GF(7741)
max_rounds = 15

println("Hydra FGLM")
println("Field: ", K)
println("Maximal number of rounds: ", max_rounds)
println("-----------------------------------------------------------")

# Trial run for precompilation
hydra_trial = Hydra_constructor(field=K, rounds_head=2)
gb_trial = extract_Hydra_gb(hydra_trial)
P_trial = parent(gb_trial[1])
bench_fglm(gb_trial, P_trial; print_time=false)
gb_trial = random_quadratic_DRL_gb(P_trial; nr_thrds=nr_thrds, print_time=false)

for rounds_head in 3:max_rounds
    println("Rounds: ", rounds_head)
    println("Hydra")
    hydra = Hydra_constructor(field=K, rounds_head=rounds_head)
    gb = extract_Hydra_gb(hydra)
    P = parent(gb[1])
    bench_fglm(gb, P)
    println("Random Quadratic")
    gb = random_quadratic_DRL_gb(P; nr_thrds=nr_thrds)
    bench_fglm(gb, P)
    println("-----------------------------------------------------------")
end
