#include "primer/hyperloglog.h"

namespace bustub {

template <typename KeyType>
HyperLogLog<KeyType>::HyperLogLog(int16_t n_bits) : cardinality_(0) {
  this->n_bits_ = n_bits;
  this->registers_ = std::vector<int16_t>(1 << n_bits_, 0);
  this->cardinality_ = 0;
}

template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeBinary(const hash_t &hash) const -> std::bitset<BITSET_CAPACITY> {
  std::bitset<BITSET_CAPACITY>bset;
  for(int i = 0; i < BITSET_CAPACITY; ++i) {
    bset[i] = static_cast<bool>((hash >> i) & 1ULL);
  }
  return bset;
}

template <typename KeyType>
auto HyperLogLog<KeyType>::PositionOfLeftmostOne(const std::bitset<BITSET_CAPACITY> &bset) const -> uint64_t {
    for (int i = n_bits_; i < BITSET_CAPACITY; ++i) {
        if (bset[BITSET_CAPACITY - 1 - i]) {  
            return i - n_bits_ + 1;
        }
    }
    return BITSET_CAPACITY - n_bits_;
}

template <typename KeyType>
auto HyperLogLog<KeyType>::AddElem(KeyType val) -> void {
  hash_t hashedVal = CalculateHash(val);
  std::bitset binaryVal = ComputeBinary(hashedVal);
  uint64_t P = PositionOfLeftmostOne(binaryVal);
  uint64_t idx = (hashedVal >> (BITSET_CAPACITY - n_bits_)) & ((1ULL << n_bits_) - 1);
  registers_[idx] = std::max(registers_[idx], static_cast<int16_t>(P));
}

template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeCardinality() -> void {
    size_t m = 1 << n_bits_;
    double sum = 0.0;
    for (size_t i = 0; i < registers_.size(); ++i) {
        sum += std::pow(2.0, -static_cast<double>(registers_[i]));
    }
    double raw_card = CONSTANT * m * m / sum; 
    cardinality_ = static_cast<size_t>(raw_card);
}


template class HyperLogLog<int64_t>;
template class HyperLogLog<std::string>;

}  
