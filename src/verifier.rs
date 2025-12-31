use crate::sis_relation::SISInstance;
use crate::sis_relation::to_principal_relation;
use labrador::iopattern::LabradorIOPattern;
use labrador::verifier::verify_principal_relation_oneround;
use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::decomposition::DecompositionFriendlySignedRepresentative;
use lattirust_arithmetic::nimue::arthur::SerArthur;
use lattirust_arithmetic::nimue::iopattern::SerIOPattern;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
use lattirust_arithmetic::traits::FromRandomBytes;
use nimue::IOPattern;
use nimue::ProofError;
use relations::Relation;
use relations::principal_relation::PrincipalRelation;
use relations::principal_relation::Witness;

pub fn verify<R: PolyRing>(proofs: &Vec<Vec<u8>>, sis: &SISInstance<R>) -> Result<(), ProofError>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
    <<R as PolyRing>::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
{
    // sis to principal relation
    println!("begin sis to principal relation");
    let (crs, index, instance, witness) = to_principal_relation(&sis);

    // verify proof
    let mut io;
    let mut arthur;
    let mut crs_use = &crs; // re-inital crs_use
    let mut index_curr = index.clone();
    let mut instance_curr = instance.clone();
    let mut round = 0;
    while crs_use.next_crs.is_some() {
        io = IOPattern::new(Box::new(format!("SIS Reduction Proof round {}", round + 1)).leak())
            .labrador_io(crs_use);
        arthur = io.to_arthur(&proofs[round]);
        (index_curr, instance_curr) =
            verify_principal_relation_oneround(&mut arthur, crs_use, &index_curr, &instance_curr)
                .expect("Failed to verify proof");
        crs_use = crs_use.next_crs.as_ref().unwrap();
        round += 1;
    }

    // final check
    io = IOPattern::new("final check").absorb_vectors::<R>(crs_use.n, crs_use.r, "witness");
    arthur = io.to_arthur(&proofs[round]);
    let s = arthur
        .next_vectors(index_curr.n, index_curr.r)
        .expect("error extracting witness");
    let witness = Witness::<R>::new(s);
    match PrincipalRelation::<R>::is_satisfied_err(&index_curr, &instance_curr, &witness) {
        Ok(_) => {
            println!("└ Verifier::verify_principal_relation: OK");
            Ok(())
        }
        Err(e) => {
            println!("└ Verifier::verify_principal_relation: ERROR {}", e);
            Err(ProofError::InvalidProof)
        }
    }
}
