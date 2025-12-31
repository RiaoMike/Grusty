#[cfg(test)]
mod tests {
    use crate::prover::prove;
    use crate::sis_relation::SISInstance;
    use crate::verifier::verify;
    use ark_std::rand;
    use ark_std::rand::Rng;
    use lattirust_arithmetic::linear_algebra::{Matrix, Vector};
    use lattirust_arithmetic::ring::PolyRing;
    use lattirust_arithmetic::ring::Pow2CyclotomicPolyRingNTT;
    use lattirust_arithmetic::ring::Zq2;
    use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;

    // Q = 2^64+1
    const Q1: u64 = 274177;
    const Q2: u64 = 67280421310721;
    pub type Z64 = Zq2<Q1, Q2>;
    const D: usize = 64;

    type R = Pow2CyclotomicPolyRingNTT<Z64, D>;

    fn random_sis_instance<R: PolyRing>(m: usize, n: usize) -> SISInstance<R>
    where
        R::BaseRing: WithSignedRepresentative,
    {
        let mut rng = rand::thread_rng();
        // this is just an example
        // let linf_bound = 3;
        // let norm_bound_squared = (linf_bound * linf_bound * n as i128 * D as i128) as f64;
        let norm_bound_squared = (Q1 as f64 * Q2 as f64).sqrt();
        let linf_bound =
            (norm_bound_squared.sqrt() / f64::sqrt((n * R::dimension()) as f64)).floor() as i128;

        let a = Matrix::rand(m, n, &mut rng);
        let s = Vector::<R>::from(
                (0..n)
                    .map(|_| {
                        R::try_from_coefficients(
                            (0..R::dimension())
                                .map(|_| {
                                    let rand_i128 = rng.sample(rand::distributions::Uniform::<i128>::new(
                                        -linf_bound,
                                        linf_bound,
                                    ));
                                    <R::BaseRing as WithSignedRepresentative>::SignedRepresentative::try_from(
                                        rand_i128
                                    ).unwrap()
                                    .into()
                                })
                                .collect::<Vec<R::BaseRing>>()
                                .as_slice(),
                        )
                        .unwrap()
                    })
                    .collect::<Vec<R>>(),
            );

        // Compute t = A * s
        let t = &a * &s;

        SISInstance::new(m, n, a, t, norm_bound_squared, s).expect("Failed to create SIS instance")
    }

    #[test]
    fn test_sis_reduction_proof() {
        // instance
        let m = 4;
        let n = 128;
        println!("begin generate random sis instance");
        let sis = random_sis_instance::<R>(m, n);
        let proofs = prove::<R>(&sis).expect("Fail to Prove");

        // print proof size
        let mut sum = 0;
        sum = proofs.iter().map(|x| x.len()).sum();
        println!(
            "Proof generated successfully! Proof size: {:.1} KB",
            sum as f64 / 1024.0
        );

        // verify proof
        let result = verify::<R>(&proofs, &sis);
        assert!(result.is_ok());
        println!("Proof verified successfully!");
    }
}
