#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(unused)]
extern crate bench_utils;
extern crate derivative;

pub mod algebra;
pub mod alt_bn128;
pub mod hashing;
pub mod merkle_tree;

pub type Error = Box<dyn std::error::Error>;