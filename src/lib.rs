use std::str::FromStr;
use uint::construct_uint;
use wasm_bindgen::prelude::*;

construct_uint! {
    struct U512(8);
}

#[derive(Copy, Clone)]
struct JacobianPoint {
    x: U512,
    y: U512,
    z: U512,
}

struct Point {
    x: U512,
    y: U512,
}

#[wasm_bindgen]
pub struct Key {
    private: U512,
    public: Point,
}

const G: Point = Point {
    x: U512([
        0x343280d6115c1d21,
        0x4a03c1d356c21122,
        0x6bb4bf7f321390b9,
        0xb70e0cbd,
        0x0,
        0x0,
        0x0,
        0x0,
    ]),
    y: U512([
        0x44d5819985007e34,
        0xcd4375a05a074764,
        0xb5f723fb4c22dfe6,
        0xbd376388,
        0x0,
        0x0,
        0x0,
        0x0,
    ]),
};
const P: U512 = U512([
    0x0000000000000001,
    0xffffffff00000000,
    0xffffffffffffffff,
    0xffffffff,
    0x0,
    0x0,
    0x0,
    0x0,
]);
const A: U512 = U512([
    0xfffffffffffffffe,
    0xfffffffeffffffff,
    0xffffffffffffffff,
    0xffffffff,
    0x0,
    0x0,
    0x0,
    0x0,
]);
const N: U512 = U512([
    0x13dd29455c5c2a3d,
    0xffff16a2e0b8f03e,
    0xffffffffffffffff,
    0xffffffff,
    0x0,
    0x0,
    0x0,
    0x0,
]);

const KEY_SIZE: usize = 224;
const WINDOW_SIZE: usize = 5;

impl U512 {
    fn inv_mod(self, m: U512) -> U512 {
        let (mut r, mut new_r) = (self % m, m);
        if r.is_zero() {
            panic!("a is not invertible");
        }
        let (mut t, mut new_t) = (U512::one(), U512::zero());
        while !new_r.is_zero() {
            let quotient = r / new_r;
            (r, new_r) = (new_r, r - quotient * new_r);
            (t, new_t) = (new_t, t + (m - quotient % m) * new_t);
            new_t %= m;
        }
        t
    }
}

impl JacobianPoint {
    fn to_affine(&self) -> Point {
        let z_inv = self.z.inv_mod(P);
        let z_inv2 = (z_inv * z_inv) % P;
        let z_inv3 = (z_inv2 * z_inv) % P;
        Point {
            x: (self.x * z_inv2) % P,
            y: (self.y * z_inv3) % P,
        }
    }

    fn zero() -> JacobianPoint {
        JacobianPoint {
            x: U512::zero(),
            y: U512::zero(),
            z: U512::one(),
        }
    }

    fn negate(&self) -> JacobianPoint {
        JacobianPoint {
            x: self.x,
            y: P - self.y,
            z: self.z,
        }
    }

    fn add(&self, other: &JacobianPoint) -> JacobianPoint {
        if self.x.is_zero() && self.y.is_zero() {
            return *other;
        }
        if other.x.is_zero() && other.y.is_zero() {
            return *self;
        }
        // https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-1998-cmo-2
        let z1z1 = (self.z * self.z) % P;
        let z2z2 = (other.z * other.z) % P;
        let u1 = (self.x * z2z2) % P;
        let u2 = (other.x * z1z1) % P;
        let s1 = (((self.y * other.z) % P) * z2z2) % P;
        let s2 = (((other.y * self.z) % P) * z1z1) % P;
        if u1 == u2 && s1 == s2 {
            return self.double();
        }
        if u1 == u2 && s1 != s2 {
            return JacobianPoint::zero();
        }
        let h = (u2 + (P - u1)) % P;
        let hh = (h * h) % P;
        let hhh = (hh * h) % P;
        let r = (s2 + (P - s1)) % P;
        let v = (u1 * hh) % P;
        let x3 = ((r * r) + (P - hhh) + (P - v) * 2) % P;
        let y3 = (r * (v + (P - x3)) + (P - s1) * hhh) % P;
        let z3 = (h * ((self.z * other.z) % P)) % P;
        JacobianPoint {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    fn double(&self) -> JacobianPoint {
        if self.x.is_zero() && self.y.is_zero() {
            return JacobianPoint::zero();
        }
        // https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-2007-bl
        let xx = (self.x * self.x) % P;
        let yy = (self.y * self.y) % P;
        let yyyy = (yy * yy) % P;
        let zz = (self.z * self.z) % P;
        let tmp_s = (self.x + yy) % P;
        let s = (((tmp_s * tmp_s) % P + (P - xx) + (P - yyyy)) * 2) % P;
        let m = (xx * 3 + A * ((zz * zz) % P)) % P;
        let t = ((m * m) + (P - s) * 2) % P;
        let x3 = t % P;
        let y3 = (m * (s + (P - t)) + (P - yyyy) * 8) % P;
        let tmp_z3 = self.z + self.y;
        let z3 = ((tmp_z3 * tmp_z3) + (P - yy) + (P - zz)) % P;
        JacobianPoint {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    fn mul(&self, other: U512) -> JacobianPoint {
        let mut result = JacobianPoint::zero();
        let mut multiplier = other % N;
        let addend = self;

        let mut precomp: [JacobianPoint; (1 << (WINDOW_SIZE - 1))] =
            [JacobianPoint::zero(); (1 << (WINDOW_SIZE - 1))];
        for i in 1..(1 << (WINDOW_SIZE - 1)) {
            precomp[i] = precomp[i - 1].add(addend);
        }

        let mut naf: [i32; KEY_SIZE + 1] = [0; KEY_SIZE + 1];
        let mask = (U512::one() << WINDOW_SIZE) - 1;
        for i in 0..naf.len() {
            if !(multiplier & U512::one()).is_zero() {
                naf[i] = (multiplier & mask).as_u32() as i32;
                if naf[i] >= (1 << (WINDOW_SIZE - 1)) {
                    naf[i] -= 1 << WINDOW_SIZE;
                    multiplier = multiplier + (-naf[i]);
                } else {
                    multiplier = multiplier - naf[i];
                }
            };
            multiplier >>= 1;
            if multiplier.is_zero() {
                break;
            }
        }
        for i in (0..KEY_SIZE + 1).rev() {
            result = result.double();
            if naf[i] == 0 {
                continue;
            }
            let ind = naf[i].unsigned_abs() as usize;
            if naf[i] > 0 {
                result = result.add(&precomp[ind]);
            } else {
                result = result.add(&precomp[ind].negate());
            }
        }

        result
    }
}

impl Point {
    fn to_jacobian(&self) -> JacobianPoint {
        JacobianPoint {
            x: self.x,
            y: self.y,
            z: U512::one(),
        }
    }
    fn mul(&self, other: U512) -> Point {
        self.to_jacobian().mul(other).to_affine()
    }
}

#[wasm_bindgen]
impl Key {
    pub fn from_private_hex(s: &str) -> Key {
        let k = U512::from_str(s).unwrap();
        let p = G.mul(k);
        Key {
            private: k,
            public: p,
        }
    }
    pub fn from_public_hex(s: &str) -> Key {
        if s.len() == (KEY_SIZE / 2 + 1) {
            let xy = U512::from_str(s).unwrap() ^ (U512::one() << (2 * KEY_SIZE + 2));
            let x = xy >> KEY_SIZE;
            let y = xy & ((U512::one() << KEY_SIZE) - 1);
            Key {
                private: U512::zero(),
                public: Point { x, y },
            }
        } else if s.len() == (KEY_SIZE / 4 + 1) {
            panic!("compressed keys not supported yet");
        } else {
            panic!("Invalid public key length");
        }
    }
    pub fn get_private_hex(&self) -> String {
        format!("{:0>64x}", self.private)
    }
    pub fn get_public_hex(&self) -> String {
        format!(
            "{:0>64x}",
            (U512::one() << (2 * KEY_SIZE + 2)) | (self.public.x << KEY_SIZE) | self.public.y
        )
    }
    pub fn derive(&self, other: &Key) -> Key {
        if !self.private.is_zero() {
            Key {
                private: U512::zero(),
                public: other.public.mul(self.private),
            }
        } else if !other.private.is_zero() {
            return Key {
                private: U512::zero(),
                public: self.public.mul(other.private),
            };
        } else {
            panic!("Cannot derive from two public keys");
        }
    }
}
