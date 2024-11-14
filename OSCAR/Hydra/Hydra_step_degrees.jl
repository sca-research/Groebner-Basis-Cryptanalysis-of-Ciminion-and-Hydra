using Oscar
include("Hydra.jl")
include("Hydra_polynomial_model.jl")


nr_thrds = 48
K = GF(7741)
max_rounds = 15
with_gb = true

println("F4 Step Degrees For Hydra")
println("Field: ", K)
println("Maximal number of rounds: ", max_rounds)
println("With Gröbner Basis: ", with_gb)
println("-----------------------------------------------------------")

for rounds_head in 3:max_rounds
    println("Rounds: ", rounds_head)
    m = 2
    hydra = Hydra_constructor(field=K, rounds_head=rounds_head)
    polys = generate_Hydra_polynomials_m_samples(hydra=hydra, m=m);
    polys = transform_Hydra_polynomial_system(hydra, polys, m);

    if with_gb
        affine_polys, polys_subs, polys_downsized_subs = non_linear_variable_substitution_Hydra_polynomial_system(hydra, polys, m; transformed=true);

        P = parent(polys_downsized_subs[1])
        variables_subs = map(i -> "x_subs_i" * string(i), 1:2 * hydra.rounds_head - 2)
        P_subs, variables_subs = polynomial_ring(hydra.field, variables_subs);
        induce(variables_subs, degrevlex(variables_subs));
        zero_vec = zero_matrix(P_subs, length(gens(P)) - length(variables_subs), 1);
        image = [vec(zero_vec[:, 1]); variables_subs];
        phi = hom(P, P_subs, image);
    
        polys_downsized_subs = map(phi, polys_downsized_subs);
    
        gb = groebner_basis_f4(ideal(polys_downsized_subs), nr_thrds=nr_thrds, info_level=2);
    else
        lin_polys = filter(poly -> total_degree(poly) == 1, polys);
        lin_polys = gens(groebner_basis_f4(ideal(lin_polys)));

        polys = filter(poly -> total_degree(poly) > 1, polys);
        polys = map(poly -> divrem(poly, lin_polys)[2], polys);

        gb = groebner_basis_f4(ideal(polys), nr_thrds=nr_thrds, info_level=2);
    end

    println("-----------------------------------------------------------")
end
