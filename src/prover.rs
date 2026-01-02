use crate::sis_relation::SISInstance;
use crate::sis_relation::to_principal_relation;
use labrador::iopattern::LabradorIOPattern;
use labrador::prover::prove_principal_relation_oneround;
use lattirust_arithmetic::challenge_set::labrador_challenge_set::LabradorChallengeSet;
use lattirust_arithmetic::challenge_set::weighted_ternary::WeightedTernaryChallengeSet;
use lattirust_arithmetic::decomposition::DecompositionFriendlySignedRepresentative;
use lattirust_arithmetic::nimue::iopattern::SerIOPattern;
use lattirust_arithmetic::nimue::merlin::SerMerlin;
use lattirust_arithmetic::ring::PolyRing;
use lattirust_arithmetic::ring::representatives::WithSignedRepresentative;
use lattirust_arithmetic::traits::FromRandomBytes;
use nimue::IOPattern;
use nimue::ProofResult;
use std::fmt::Debug;

pub fn prove<R: PolyRing>(sis: &SISInstance<R>) -> ProofResult<Vec<Vec<u8>>>
where
    LabradorChallengeSet<R>: FromRandomBytes<R>,
    WeightedTernaryChallengeSet<R>: FromRandomBytes<R>,
    <R as PolyRing>::BaseRing: WithSignedRepresentative,
    <R::BaseRing as WithSignedRepresentative>::SignedRepresentative:
        DecompositionFriendlySignedRepresentative,
    <R as TryFrom<u128>>::Error: Debug,
{
    //init
    let mut proofs: Vec<Vec<u8>> = Vec::new();
    // sis to principal relation
    println!("begin sis to principal relation");
    let (crs, index, instance, witness) = to_principal_relation(&sis);
    println!("CRS Structure Info:");
    println!("proof_size: {} KB", crs.proof_size() as f64 / 1024.0);
    println!(
        "last_prover_message_size: {} KB",
        crs.last_prover_message_size() as f64 / 1024.0
    );

    let mut crs_use = &crs;

    // begin proof
    println!("begin proof");
    let mut io;
    let mut merlin;

    let mut index_curr = index.clone();
    let mut instance_curr = instance.clone();
    let mut witness_curr = witness.clone();
    let mut round = 1;
    while crs_use.next_crs.is_some() {
        println!("round {}", round);
        io = IOPattern::new(Box::new(format!("SIS Reduction Proof round {}", round)).leak())
            .labrador_io(crs_use);
        merlin = io.to_merlin();
        (index_curr, instance_curr, witness_curr) = prove_principal_relation_oneround(
            &mut merlin,
            crs_use,
            &index_curr,
            &instance_curr,
            &witness_curr,
        )
        .expect("Failed to prove principal relation");
        proofs.push(merlin.transcript().to_vec());
        crs_use = crs_use.next_crs.as_ref().unwrap();
        round += 1;
    }
    // 最后一轮，直接添加 witness
    io = IOPattern::new("final check").absorb_vectors::<R>(crs_use.n, crs_use.r, "witness");
    merlin = io.to_merlin();
    merlin
        .absorb_vectors::<R>(&witness_curr.s)
        .expect("error absorbing witness");
    proofs.push(merlin.transcript().to_vec());
    Ok(proofs)
}
