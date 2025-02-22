/**
 * @file    karp_rabin_hashing.cpp
 * @section LICENCE
 *
 * This file is part of Lazy-AVLG v0.1.0
 * See: https://github.com/dominikkempa/lz77-to-slp
 *
 * Copyright (C) 2021
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/
 //WZ:Modified for the purpose of computing hashes of subsequent factors of length k faster

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <random>
#include "utils.h"
#include "krfp.h"

using namespace std;

namespace karp_rabin_hashing {

//=============================================================================
// Base and exponent used in Karp-Rabin hashing.
//=============================================================================
std::uint64_t hash_variable;
std::uint64_t hash_power_k; // hash variable to power k, where k is passed during initialisation
std::vector<std::uint64_t> hash_power_table; // precomputed hash variable powers up to power k for faster operations
std::uint64_t mersenne_prime_exponent;
std::uint64_t inverse;

//=============================================================================
// Return (a * b) mod p, where p = (2^k) - 1.
// Requires a, b <= 2^k. Tested for k = 1, .., 63.
//=============================================================================
inline std::uint64_t mul_mod_mersenne(
    const std::uint64_t a,
    const std::uint64_t b,
    const std::uint64_t k) {
  const std::uint64_t p = ((std::uint64_t)1 << k) - 1;
  __extension__ const unsigned __int128 ab =
    (unsigned __int128)a *
    (unsigned __int128)b;
  std::uint64_t lo = (std::uint64_t)ab;
  const std::uint64_t hi = (ab >>  (uint64_t) 64);
  lo = (lo & p) + ((lo >> k) + (hi << ( (uint64_t) 64 - k)));
  lo = (lo & p) + (lo >> k);
  return lo == p ?  (uint64_t) 0 : lo;
}

//=============================================================================
// Return a mod p, where p = (2^k) - 1.
// Works for any a in [0..2^64).
// Tested for k = 1, .., 63.
//=============================================================================
inline std::uint64_t mod_mersenne(
    std::uint64_t a,
    const std::uint64_t k) {
  std::uint64_t p = ((std::uint64_t)1 << k) -  (uint64_t) 1;
  if (k < (uint64_t) 32) {

    // We need to check if a <= 2^(2k).
    const std::uint64_t threshold = ((std::uint64_t)1 << (k <<  (uint64_t) 1));
    if (a <= threshold) {
      a = (a & p) + (a >> k);
      a = (a & p) + (a >> k);
      return a == p ?  (uint64_t) 0 : a;
    } else return a % p;
  } else {

    // We are guaranteed that a < 2^(2k)
    // because a < 2^64 <= 2^(2k).
    a = (a & p) + (a >> k);
    a = (a & p) + (a >> k);
    return a == p ?  (uint64_t) 0 : a;
  }
}

//=============================================================================
// Return random number x in [0..p), where p = (2^k) - 1.
//=============================================================================
std::uint64_t rand_mod_mersenne(const std::uint64_t k) {
  const std::uint64_t p = ((std::uint64_t)1 << k) -  (uint64_t) 1;
  return utils::random_int<std::uint64_t>(
      (std::uint64_t)0, (std::uint64_t(p -  (uint64_t) 1)));
}

//=============================================================================
// Return (a^n) mod p, where p = (2^k) - 1.
//=============================================================================
std::uint64_t  pow_mod_mersenne(
    const std::uint64_t a,
    std::uint64_t n,
    const std::uint64_t k) {
  std::uint64_t pow = mod_mersenne(a, k);
  std::uint64_t ret = mod_mersenne( (uint64_t) 1, k);
  while (n >  (uint64_t) 0) {
    if (n &  (uint64_t) 1)
      ret = mul_mod_mersenne(ret, pow, k);
    pow = mul_mod_mersenne(pow, pow, k);
    n >>=  (uint64_t) 1;
  }
  return ret;
}

//=============================================================================
// Given Karp-Rabin hashes of two substrings, return
// the Karp-Rabin hash of their concatenation.
//=============================================================================
std::uint64_t concat(
    const std::uint64_t left_hash,
    const std::uint64_t right_hash,
    const std::uint64_t right_len) {
  const std::uint64_t pow = pow_mod_mersenne(
      hash_variable, right_len, mersenne_prime_exponent);
  const std::uint64_t tmp = mul_mod_mersenne(
      left_hash, pow, mersenne_prime_exponent);
  const std::uint64_t ret = mod_mersenne(
      tmp + right_hash, mersenne_prime_exponent);
  return ret;
}

std::uint64_t subtract(
    const std::uint64_t long_hash,
    const std::uint64_t short_hash,
    const std::uint64_t right_len) {
  const std::uint64_t pow = pow_mod_mersenne(
      hash_variable, right_len, mersenne_prime_exponent);
  const std::uint64_t tmp = mul_mod_mersenne(
      short_hash, pow, mersenne_prime_exponent);
  const std::uint64_t p = ((std::uint64_t)1 << mersenne_prime_exponent) - 1;
  return (long_hash >= tmp) ?
    (long_hash - tmp) :
    ((long_hash + p) - tmp);
}


std::uint64_t leftshift(const std::uint64_t hash){
	return mul_mod_mersenne(hash, inverse, mersenne_prime_exponent);
}

std::uint64_t concat_k(// here we assume, that right_len is a constant passed during initialization
    const std::uint64_t left_hash,
    const std::uint64_t right_hash) {
  const std::uint64_t tmp = mul_mod_mersenne(
      left_hash, hash_power_k, mersenne_prime_exponent);
  const std::uint64_t ret = mod_mersenne(
      tmp + right_hash, mersenne_prime_exponent);
  return ret;
}

std::uint64_t fast_concat(// faster computation (O(1) instead of O(log (right_len))) assuming right_len <= k
    const std::uint64_t left_hash,
    const std::uint64_t right_hash,
    const std::uint64_t right_len) {
  const std::uint64_t pow = hash_power_table[right_len];
  const std::uint64_t tmp = mul_mod_mersenne(
      left_hash, pow, mersenne_prime_exponent);
  const std::uint64_t ret = mod_mersenne(
      tmp + right_hash, mersenne_prime_exponent);
  return ret;
}


std::uint64_t subtract_k( // here we assume, that right_len is a constant passed during initialization
    const std::uint64_t long_hash,
    const std::uint64_t short_hash) {
  const std::uint64_t tmp = mul_mod_mersenne(
      short_hash, hash_power_k, mersenne_prime_exponent);
  const std::uint64_t p = ((std::uint64_t)1 << mersenne_prime_exponent) - 1;
  return (long_hash >= tmp) ?
    (long_hash - tmp) :
    ((long_hash + p) - tmp);
}

std::uint64_t fast_subtract(// faster computation (O(1) instead of O(log (right_len))) assuming right_len <= k
    const std::uint64_t long_hash,
    const std::uint64_t short_hash,
    const std::uint64_t right_len) {
  const std::uint64_t pow = hash_power_table[right_len];
  const std::uint64_t tmp = mul_mod_mersenne(
      short_hash, pow, mersenne_prime_exponent);
  const std::uint64_t p = ((std::uint64_t)1 << mersenne_prime_exponent) - 1;
  return (long_hash >= tmp) ?
    (long_hash - tmp) :
    ((long_hash + p) - tmp);
}


//computation of the inverse of hash_variable
void compute_inverse(){
	inverse=1;
	std::uint64_t multip;
	std::uint64_t p = ((std::uint64_t)1 << mersenne_prime_exponent) - 1;
	std::uint64_t a = hash_variable;
	while(a!=1){
		multip = p/a;
		multip=mul_mod_mersenne(multip, p-1,mersenne_prime_exponent);
		inverse =mul_mod_mersenne(inverse,multip,mersenne_prime_exponent);
		a=p%a;
	}
}




//=============================================================================
// Initialize the base and exponent for Karp-Rabin hashing.
//=============================================================================
void init(uint64_t k) {
  hash_power_table.clear(); // clean previous initialisation
  mersenne_prime_exponent = 61; //do not change this
  //hash_variable = 838321722101173966; // fixed random value to remove variance in tests
  hash_variable = rand_mod_mersenne(mersenne_prime_exponent);
  compute_inverse();
  hash_power_table.push_back(mod_mersenne( (uint64_t) 1, k));
  for(uint64_t i = 1 ; i <= k; ++i){
	  hash_power_table.push_back(mul_mod_mersenne(
      hash_power_table[i - 1], hash_variable, mersenne_prime_exponent));
  }
  hash_power_k = pow_mod_mersenne(hash_variable, k, mersenne_prime_exponent);
}
}

