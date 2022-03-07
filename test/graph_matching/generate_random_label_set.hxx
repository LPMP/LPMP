#pragma once

#include <vector>
#include <algorithm>
#include <random>

std::vector<std::size_t> generate_random_label_set(const std::size_t no_labels, const std::size_t max_label)
{
   assert(no_labels <= max_label);
   std::vector<std::size_t> labels(max_label);
   std::iota(labels.begin(), labels.end(), 0);

   // randomly shuffle, cut off last entries, sort
   std::random_device rd;
   std::mt19937 g(rd());
   std::shuffle(labels.begin(), labels.end(), g);
   labels.resize(no_labels);
   std::sort(labels.begin(), labels.end());
   return labels;
}
