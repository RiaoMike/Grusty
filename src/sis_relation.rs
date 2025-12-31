use labrador::common_reference_string::CommonReferenceString;
use lattirust_arithmetic::linear_algebra::Matrix;
use lattirust_arithmetic::linear_algebra::SymmetricMatrix;
use lattirust_arithmetic::linear_algebra::Vector;
use lattirust_arithmetic::ring::PolyRing;
use num_traits::Zero;
use relations::principal_relation::Size;
use relations::principal_relation::{ConstantQuadraticConstraint, QuadraticConstraint};
use relations::principal_relation::{Index, Instance, Witness};

use lattirust_arithmetic::traits::WithL2Norm;
use num_bigint::BigUint;
use num_traits::cast::FromPrimitive;

#[derive(Debug)]
pub struct SISInstance<R: PolyRing> {
    // size of row of A
    pub m: usize,
    // size of col of A
    pub n: usize,
    // commitment matrix A
    pub a: Matrix<R>,
    // ..
    pub t: Vector<R>,
    // norm bound squared of L2
    pub norm_bound_squared: f64,

    // private s
    pub s: Vector<R>,
}

impl<R: PolyRing> SISInstance<R> {
    pub fn new(
        m: usize,
        n: usize,
        a: Matrix<R>,
        t: Vector<R>,
        norm_bound_squared: f64,
        s: Vector<R>,
    ) -> Result<Self, &'static str> {
        let instance: Self = Self {
            m,
            n,
            a,
            t,
            norm_bound_squared,
            s,
        };
        if !instance.is_satisfied() {
            Err("SISInstance is not satisfied")
        } else {
            Ok(instance)
        }
    }

    fn is_satisfied(&self) -> bool {
        self.m == self.t.len()
            && self.n == self.s.len()
            && self.a.nrows() == self.m
            && self.a.ncols() == self.n
            && self.s.l2_norm_squared() <= BigUint::from_f64(self.norm_bound_squared).unwrap()
            && &self.a * &self.s == self.t
    }
}

pub fn to_principal_relation<R: PolyRing>(
    instance: &SISInstance<R>,
) -> (CommonReferenceString<R>, Index<R>, Instance<R>, Witness<R>) {
    let n = instance.n;
    let target_num_witnesses = ((n as f64).powf(1.0 / 3.0).trunc() as usize).max(1);

    let num_witnesses = (1..=target_num_witnesses)
        .rev()
        .find(|&r| n % r == 0)
        .unwrap_or(1);
    let witness_len = n / num_witnesses;

    let size = Size {
        num_witnesses,
        witness_len,
        norm_bound_sq: instance.norm_bound_squared,
        num_constraints: instance.m,
        num_constant_constraints: 1,
    };
    let crs = CommonReferenceString::new_for_size(size);
    let index = Index::new(&size);
    let inst = Instance::<R> {
        quad_dot_prod_funcs: (0..index.num_constraints)
            .map(|i| {
                let row_i = instance.a.row(i);
                let split_row: Vec<Vector<R>> = row_i
                    .iter()
                    .cloned()
                    .collect::<Vec<_>>()
                    .chunks(witness_len)
                    .map(|chunk| Vector::from_vec(chunk.to_vec()))
                    .collect();
                QuadraticConstraint::new_linear(
                    // SymmetricMatrix::<R>::zero(num_witnesses),
                    split_row,
                    instance.t[i].clone(),
                )
            })
            .collect(),
        ct_quad_dot_prod_funcs: vec![ConstantQuadraticConstraint::<R>::new(
            SymmetricMatrix::<R>::zero(num_witnesses),
            vec![Vector::zeros(witness_len); num_witnesses],
            R::BaseRing::zero(),
        )],
    };

    let witness_splits: Vec<Vector<R>> = instance
        .s
        .iter()
        .cloned()
        .collect::<Vec<_>>()
        .chunks(witness_len)
        .map(|chunk| Vector::from_vec(chunk.to_vec()))
        .collect();
    let witness = Witness::new(witness_splits);

    (crs, index, inst, witness)
}
