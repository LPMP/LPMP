#ifndef LPMP_COMPUTE_RANDOM_LABEL_SET
#define LPMP_COMPUTE_RANDOM_LABEL_SET

#include <vector>

std::vector<std::size_t> generate_random_label_set(const std::size_t no_labels, const std::size_t max_label)
{
   assert(no_labels <= max_label);
   std::vector<std::size_t> labels(max_label);
   std::iota(labels.begin(), labels.end(), 0);

   // randomly shuffle, cut off last entries, sort
   std::random_shuffle(labels.begin(), labels.end());
   labels.resize(no_labels);
   std::sort(labels.begin(), labels.end());
   return labels;
}

#endif // LPMP_COMPUTE_RANDOM_LABEL_SET

