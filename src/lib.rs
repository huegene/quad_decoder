//#![doc = include_str!("../README.md")]
#![no_std]
// The quad decoder module combines a few data structures and algorithms to enable extraction of
// distance data from a quadrature current signal.

use core::ops::{Add, Mul, Sub};
use fixed::types::*;
use heapless::Vec;
#[allow(unused_macros)]
macro_rules! fixed {
    ($val:expr) => {
        I16F16::from_num($val)
    };
}

impl Vertex {
    pub fn new(x: I16F16, y: I16F16) -> Self {
        Self {
            x: fixed!(x),
            y: fixed!(y),
        }
    }
    pub fn abs(&self) -> I16F16 {
        let xx: I32F32 = I32F32::from_num(self.x) * I32F32::from_num(self.x);
        let yy: I32F32 = I32F32::from_num(self.y) * I32F32::from_num(self.y);
        let sqr = xx + yy;
        I16F16::from_num(cordic::sqrt(sqr))
    }
}

impl Circle {
    pub fn new(x: I16F16, y: I16F16, r: I16F16) -> Self {
        Self {
            x: fixed!(x),
            y: fixed!(y),
            r: fixed!(r),
        }
    }
    pub fn exp_filt(self, rhs: Circle, alpha: I16F16) -> Self {
        let xx = self.x * (fixed!(1.0) - alpha) + rhs.x * alpha;
        let yy = self.y * (fixed!(1.0) - alpha) + rhs.y * alpha;
        let rr = self.r * (fixed!(1.0) - alpha) + rhs.r * alpha;
        //self.x=xx;
        //self.y=yy;
        //self.r=rr;
        Self {
            x: xx,
            y: yy,
            r: rr,
        }
    }
}

/// Basic point structure, for now basic required traits are just derived
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Vertex {
    pub x: I16F16,
    pub y: I16F16,
}

/// Basic circle structure, same as Vertex, just with radius. Maybe some inheritance should be
/// considered
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct Circle {
    pub x: I16F16,
    pub y: I16F16,
    pub r: I16F16,
}

impl Mul<I16F16> for Vertex {
    // The multiplication of rational numbers is a closed operation.
    type Output = Self;

    fn mul(self, rhs: I16F16) -> Self {
        let xx = self.x * rhs;
        let yy = self.y * rhs;
        Vertex { x: xx, y: yy }
    }
}

impl Add for Vertex {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub for Vertex {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Add<Vertex> for Circle {
    type Output = Self;

    fn add(self, other: Vertex) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            r: self.r,
        }
    }
}

pub trait ExpFilt<T, F> {
    fn exp_filt(self, rhs: T, alpha: F) -> Self;
}

/// Simple formula to calculate center point and radius of a circle from 2 points
#[allow(dead_code)]
pub fn circle_from_three_vertex(vertex: &Vec<Vertex, 3>) -> Circle {
    // implemented like in from https://www.geeksforgeeks.org/equation-of-circle-when-three-points-on-the-circle-are-given/
    let bx1 = I64F0::from_num(vertex[0].x);
    let bx2 = I64F0::from_num(vertex[1].x);
    let bx3 = I64F0::from_num(vertex[2].x);
    let by1 = I64F0::from_num(vertex[0].y);
    let by2 = I64F0::from_num(vertex[1].y);
    let by3 = I64F0::from_num(vertex[2].y);
    let x12 = bx1 - bx2;
    let x13 = bx1 - bx3;
    let y12 = by1 - by2;
    let y13 = by1 - by3;

    let x31 = bx3 - bx1;
    let x21 = bx2 - bx1;
    let y31 = by3 - by1;
    let y21 = by2 - by1;

    let sx13 = bx1 * bx1 - bx3 * bx3;
    let sy13 = by1 * by1 - by3 * by3;
    let sx21 = bx2 * bx2 - bx1 * bx1;
    let sy21 = by2 * by2 - by1 * by1;

    let f = (sx13 * x12 + sy13 * x12 + sx21 * x13 + sy21 * x13) / (2 * (y31 * x12 - y21 * x13));

    let g = (sx13 * y12 + sy13 * y12 + sx21 * y13 + sy21 * y13) / (2 * (x31 * y12 - x21 * y13));
    let c = -bx1 * bx1 - by1 * by1 - 2 * g * bx1 - 2 * f * by1;
    let sqr_of_r = g * g + f * f - c;

    Circle {
        x: I16F16::from_num(-g),
        y: I16F16::from_num(-f),
        r: I16F16::from_num(cordic::sqrt(sqr_of_r)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_circle() {
        let mut buffer = Vec::<Vertex, 3>::new();
        //println!("circle data:");
        let _ = buffer.push(Vertex {
            x: fixed!(1.0),
            y: fixed!(0.0),
        });
        let _ = buffer.push(Vertex {
            x: fixed!(0.0),
            y: fixed!(1.0),
        });
        let _ = buffer.push(Vertex {
            x: fixed!(-1.0),
            y: fixed!(0.0),
        });
        let circ = circle_from_three_vertex(&buffer);
        assert_eq!(
            circ,
            Circle {
                x: fixed!(0.0),
                y: fixed!(0.0),
                r: fixed!(1.0)
            }
        );
        //println!("{:?}",circ);
    }

    #[test]
    fn test_circle_with_offset_and_radius() {
        let mut buffer = Vec::<Vertex, 3>::new();
        //println!("circle data:");
        let _ = buffer.push(Vertex {
            x: fixed!(5.0),
            y: fixed!(2.0),
        });
        let _ = buffer.push(Vertex {
            x: fixed!(2.0),
            y: fixed!(5.0),
        });
        let _ = buffer.push(Vertex {
            x: fixed!(-1.0),
            y: fixed!(2.0),
        });
        let circ = circle_from_three_vertex(&buffer);
        assert_eq!(
            circ,
            Circle {
                x: fixed!(2.0),
                y: fixed!(2.0),
                r: fixed!(3.0)
            }
        );
        //println!("{:?}",circ);
    }
    #[test]
    fn test_circle_with_offset() {
        let mut buffer = Vec::<Vertex, 3>::new();
        //println!("circle data:");
        let _ = buffer.push(Vertex {
            x: fixed!(3.0),
            y: fixed!(2.0),
        });
        let _ = buffer.push(Vertex {
            x: fixed!(2.0),
            y: fixed!(3.0),
        });
        let _ = buffer.push(Vertex {
            x: fixed!(1.0),
            y: fixed!(2.0),
        });
        let circ = circle_from_three_vertex(&buffer);
        assert_eq!(
            circ,
            Circle {
                x: fixed!(2.0),
                y: fixed!(2.0),
                r: fixed!(1.0)
            }
        );
        //println!("{:?}",circ);
    }
    #[test]
    fn test_scalar_product() {
        let a = Vertex {
            x: fixed!(1.0),
            y: fixed!(1.0),
        };
        let c = fixed!(2.0);
        //println!("{:?}", a*c);
        let b = Vertex {
            x: fixed!(2.0),
            y: fixed!(2.0),
        };
        assert_eq!(b, a * c);
    }
    #[test]
    fn test_addition_for_vertex() {
        let a = Vertex {
            x: fixed!(1.0),
            y: fixed!(2.0),
        };
        let b = Vertex {
            x: fixed!(3.0),
            y: fixed!(4.0),
        };
        let c = Vertex {
            x: fixed!(4.0),
            y: fixed!(6.0),
        };
        //println!("{:?}", a*c);
        assert_eq!(c, a + b);
    }

    #[test]
    fn test_exp_filter_for_circle() {
        let circ_old = Circle {
            x: fixed!(0.0),
            y: fixed!(0.0),
            r: fixed!(0.0),
        };
        let circ_new = Circle {
            x: fixed!(2000.0),
            y: fixed!(2000.0),
            r: fixed!(2000.0),
        };
        let alpha = fixed!(0.5);
        let circ_filt1 = Circle {
            x: fixed!(1000.0),
            y: fixed!(1000.0),
            r: fixed!(1000.0),
        };
        let circ_filt3 = Circle {
            x: fixed!(1500),
            y: fixed!(1500),
            r: fixed!(1500),
        };
        let circ_filt2 = circ_old.exp_filt(circ_new, alpha);
        assert_eq!(circ_filt1, circ_filt2);
        let circ_filt2 = circ_filt2.exp_filt(circ_new, alpha);
        assert_eq!(circ_filt3, circ_filt2);
    }

    #[test]
    fn test_addition_for_circle_and_vertex() {
        let a = Vertex {
            x: fixed!(1.0),
            y: fixed!(2.0),
        };
        let b = Circle {
            x: fixed!(3.0),
            y: fixed!(4.0),
            r: fixed!(1.0),
        };
        let c = Circle {
            x: fixed!(4.0),
            y: fixed!(6.0),
            r: fixed!(1.0),
        };
        //println!("{:?}", a*c);
        assert_eq!(c, b + a);
    }
}
