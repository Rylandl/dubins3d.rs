use core::f64::consts::PI as PI;
use core::cmp::Ordering;
use core::fmt;

use crate::mod2pi;

#[derive(Copy, Clone, PartialEq, Eq)]
pub enum SegmentType {
    RIGHT,
    STRAIGHT,
    LEFT,
    NONE
}

impl fmt::Debug for SegmentType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SegmentType::RIGHT => write!(f, "R"),
            SegmentType::STRAIGHT => write!(f, "S"),
            SegmentType::LEFT => write!(f, "L"),
            SegmentType::NONE => write!(f, "X"),
       }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct ManeuverCase {
    pub a: SegmentType,
    pub b: SegmentType,
    pub c: SegmentType
}

#[derive(Debug, Copy, Clone)]
pub struct DubinsStruct {
    pub t: f64,
    pub p: f64,
    pub q: f64,
    pub length: f64,
    pub case: ManeuverCase
}

pub struct DubinsManeuver2D {
    pub qi: (f64, f64, f64),
    pub qf: (f64, f64, f64),
    pub rhomin: f64,
    pub maneuver: DubinsStruct
}

// pub fn mod2pi(th: f64) -> f64{
//     let t = th % (2.0 * PI);
//     if t < 0.0 {
//         return t + 2.0 * PI;
//     }
//     return t;
// }

impl DubinsManeuver2D {
    pub fn new(qi: (f64, f64, f64), qf: (f64, f64, f64), rhomin: f64, min_length: f64, disable_ccc: bool) -> Self {
        let mut maneuver = DubinsManeuver2D{
            qi: qi, 
            qf: qf,
            rhomin: rhomin,
            maneuver: DubinsStruct{
                        t: 0.0,
                        p: 0.0,
                        q: 0.0, 
                        length: core::f64::INFINITY,
                        case: ManeuverCase{a: SegmentType::NONE, b: SegmentType::NONE, c: SegmentType::NONE}
                    }
        };
    
        let dx = maneuver.qf.0 - maneuver.qi.0;
        let dy = maneuver.qf.1 - maneuver.qi.1;
        let d = (dx*dx + dy*dy).sqrt() / maneuver.rhomin;
    
        // Normalize the problem using rotation
        let rotation_angle = mod2pi(dy.atan2(dx));
        let a = mod2pi(maneuver.qi.2 - rotation_angle);
        let b = mod2pi(maneuver.qf.2 - rotation_angle);
    
        let (sa, ca) = (a.sin(),  a.cos());
        let (sb, cb) = (b.sin(),  b.cos());
    
        let path_lsl = _lsl(&maneuver, a, b, d, sa, ca, sb, cb);
        let path_rsr = _rsr(&maneuver, a, b, d, sa, ca, sb, cb);
        let path_lsr = _lsr(&maneuver, a, b, d, sa, ca, sb, cb);
        let path_rsl = _rsl(&maneuver, a, b, d, sa, ca, sb, cb);
        
        let path_rlr = _rlr(&maneuver, a, b, d, sa, ca, sb, cb);
        let path_lrl = _lrl(&maneuver, a, b, d, sa, ca, sb, cb);
        
        let mut _paths: Vec<DubinsStruct> = Vec::new();
        if disable_ccc {
            _paths.extend([path_lsl, path_rsr, path_lsr, path_rsl]);
        }
        else {
            _paths.extend([path_lsl, path_rsr, path_lsr, path_rsl, path_rlr, path_lrl]);
        }
    
        let thresh = maneuver.rhomin * 1e-5;
        if d.abs() < thresh && a.abs() < thresh && b.abs() < thresh {
            let dist_2d = (maneuver.qi.0 - maneuver.qf.0).abs()
                            .max((maneuver.qi.1 - maneuver.qf.1).abs());
            if dist_2d < thresh {
                let path_c = _c(&maneuver);
                _paths.clear();
                _paths.extend([path_c]);
            }
        }
        
        _paths.sort_by(|p1, p2| {
            let a = p1.length;
            let b = p2.length;
            match (a.is_nan(), b.is_nan()) {
                (true, true) => Ordering::Equal,
                (true, false) => Ordering::Greater,
                (false, true) => Ordering::Less,
                (false, false) => a.partial_cmp(&b).unwrap(),
            }
        });
       
        if min_length == core::f64::NEG_INFINITY {
            maneuver.maneuver = _paths[0];
        }
        else {
            for path in _paths {
                if path.length >= min_length {
                    maneuver.maneuver = path;
                    break;
                }
            }
        }
        return maneuver;
    }
}   

fn _lsl(maneuver: &DubinsManeuver2D, a: f64, b: f64, d: f64, sa: f64, ca: f64, sb: f64, cb: f64) -> DubinsStruct {
    let aux = (cb - ca).atan2(d + sa - sb);
    let t  = mod2pi(-a + aux);
    let p = (2.0 + d*d - 2.0 * (a-b).cos() + 2.0 * d * (sa - sb)).sqrt();
    let q = mod2pi(b - aux);
    let length = (t+p+q) * maneuver.rhomin;
    
    let ds = DubinsStruct {
        t: t,
        p: p,
        q: q,
        length: length,
        case: ManeuverCase {a: SegmentType::LEFT, b: SegmentType::STRAIGHT, c: SegmentType::LEFT}
    };
    
    return ds;
}

fn _rsr(maneuver: &DubinsManeuver2D, a: f64, b: f64, d: f64, sa: f64, ca: f64, sb: f64, cb: f64) -> DubinsStruct {
    let aux = (ca - cb).atan2(d - sa + sb);
    let t  = mod2pi(a - aux);
    let p = (2.0 + d*d - 2.0 * (a-b).cos() + 2.0 * d * (sb - sa)).sqrt();
    let q = mod2pi(mod2pi(-b) + aux);
    let length = (t+p+q) * maneuver.rhomin;
    
    let ds = DubinsStruct {
        t: t,
        p: p,
        q: q,
        length: length,
        case: ManeuverCase {a: SegmentType::RIGHT, b: SegmentType::STRAIGHT, c: SegmentType::RIGHT}
    };
    
    return ds;
}

fn _lsr(maneuver: &DubinsManeuver2D, a: f64, b: f64, d: f64, sa: f64, ca: f64, sb: f64, cb: f64) -> DubinsStruct {
    let aux1 = -2.0 + d*d + 2.0 * (a-b).cos() + 2.0 * d * (sa + sb);
    let mut t = core::f64::INFINITY;
    let mut q = core::f64::INFINITY;
    let mut p = core::f64::INFINITY;
    if aux1 > 0.0 {
        p = aux1.sqrt();
        let aux2 = (-ca-cb).atan2(d+sa+sb) - (-2.0/p).atan();
        t = mod2pi(-a + aux2);
        q = mod2pi(-mod2pi(b) + aux2);
    }

    let length = (t+p+q) * maneuver.rhomin;

    let ds = DubinsStruct {
        t: t,
        p: p,
        q: q,
        length: length,
        case: ManeuverCase {a: SegmentType::LEFT, b: SegmentType::STRAIGHT, c: SegmentType::RIGHT}
    };
    
    return ds;
}

fn _rsl(maneuver: &DubinsManeuver2D, a: f64, b: f64, d: f64, sa: f64, ca: f64, sb: f64, cb: f64) -> DubinsStruct {
    let aux1 = d*d - 2.0 + 2.0 * (a-b).cos() - 2.0 * d * (sa + sb);
    let mut t = core::f64::INFINITY;
    let mut q = core::f64::INFINITY;
    let mut p = core::f64::INFINITY;
    if aux1 > 0.0 {
        p = aux1.sqrt();
        let aux2 = (ca+cb).atan2(d-sa-sb) - (2.0/p).atan();
        t = mod2pi(a - aux2);
        q = mod2pi(mod2pi(b) - aux2);
    }

    let length = (t+p+q) * maneuver.rhomin;

    let ds = DubinsStruct {
        t: t,
        p: p,
        q: q,
        length: length,
        case: ManeuverCase {a: SegmentType::RIGHT, b: SegmentType::STRAIGHT, c: SegmentType::LEFT}
    };
    
    return ds;
}

fn _rlr(maneuver: &DubinsManeuver2D, a: f64, b: f64, d: f64, sa: f64, ca: f64, sb: f64, cb: f64) -> DubinsStruct {
    let aux = (6.0 - d*d + 2.0 * (a-b).cos() + 2.0 * d * (sa-sb))/8.0;
    let mut t = core::f64::INFINITY;
    let mut q = core::f64::INFINITY;
    let mut p = core::f64::INFINITY;
    if aux.abs() <= 1.0 {
        p = mod2pi(-aux.acos());
        t = mod2pi(a - (ca-cb).atan2(d-sa+sb) + p/2.0);
        q = mod2pi(a - b - t + p);
    }

    let length = (t+p+q) * maneuver.rhomin;
    let ds = DubinsStruct {
        t: t,
        p: p,
        q: q,
        length: length,
        case: ManeuverCase {a: SegmentType::RIGHT, b: SegmentType::LEFT, c: SegmentType::RIGHT}
    };
    
    return ds;
}

fn _lrl(maneuver: &DubinsManeuver2D, a: f64, b: f64, d: f64, sa: f64, ca: f64, sb: f64, cb: f64) -> DubinsStruct {
    let aux = (6.0 - d*d + 2.0 * (a-b).cos() + 2.0 * d * (-sa+sb))/8.0;
    let mut t = core::f64::INFINITY;
    let mut q = core::f64::INFINITY;
    let mut p = core::f64::INFINITY;
    if aux.abs() <= 1.0 {
        p = mod2pi(-aux.acos());
        t = mod2pi(-a + (-ca+cb).atan2(d+sa-sb) + p/2.0);
        q = mod2pi(b - a - t + p);
    }

    let length = (t+p+q) * maneuver.rhomin;
    let ds = DubinsStruct {
        t: t,
        p: p,
        q: q,
        length: length,
        case: ManeuverCase {a: SegmentType::LEFT, b: SegmentType::RIGHT, c: SegmentType::LEFT}
    };
    
    return ds;
}

fn _c(maneuver: &DubinsManeuver2D) -> DubinsStruct{
    let ds = DubinsStruct {
        t: 0.0,
        p: 2.0*PI,
        q: 0.0,
        length: 2.0*PI * maneuver.rhomin,
        case: ManeuverCase {a: SegmentType::RIGHT, b: SegmentType::RIGHT, c: SegmentType::RIGHT}
    };
    
    return ds;
}

pub fn get_coordinates_at(maneuver: &DubinsManeuver2D, offset: f64) -> (f64,f64,f64) {
    let n_offset = offset / maneuver.rhomin;

    let qi: (f64,f64,f64) = (0.0, 0.0, maneuver.qi.2);

    let l1 = maneuver.maneuver.t;
    let l2 = maneuver.maneuver.p;
    let q1 = get_position_in_segment(l1, qi, maneuver.maneuver.case.a);
    let q2 = get_position_in_segment(l2, q1, maneuver.maneuver.case.b);

    let mut q: (f64, f64, f64);
    if n_offset < l1 {
        q = get_position_in_segment(n_offset, qi, maneuver.maneuver.case.a);
    }
    else if n_offset < (l1 + l2) {
        q = get_position_in_segment(n_offset - l1, q1, maneuver.maneuver.case.b);
    }
    else {
        q = get_position_in_segment(n_offset - l1 - l2, q2, maneuver.maneuver.case.c);
    }
    q.0 = q.0 * maneuver.rhomin + maneuver.qi.0;
    q.1 = q.1 * maneuver.rhomin + maneuver.qi.1;
    q.2 = mod2pi(q.2);

    return q;
}

fn get_position_in_segment(offset: f64, qi: (f64,f64,f64), case: SegmentType) -> (f64,f64,f64){
    let mut q: (f64, f64, f64) = (0.0, 0.0, 0.0);
    if case == SegmentType::LEFT {
        q = (
            qi.0 + (qi.2 + offset).sin() - qi.2.sin(),
            qi.1 - (qi.2 + offset).cos() + qi.2.cos(),
            qi.2 + offset
        );
    }
    else if case == SegmentType::RIGHT {
        q = (
            qi.0 - (qi.2 - offset).sin() + qi.2.sin(),
            qi.1 + (qi.2 - offset).cos() - qi.2.cos(),
            qi.2 - offset
        );
    }
    else if case == SegmentType::STRAIGHT {
        q = (
            qi.0 + (qi.2).cos() * offset,
            qi.1 + (qi.2).sin() * offset,
            qi.2
        );
    }
    return q;
}