seliNDRIx: Your Package Description Here
================

buld+\_b— output: github_document —

<!-- README.md is generated from README.Rmd. Please edit that file -->

# seliNDRIx

<!-- badges: start -->
<!-- badges: end -->

The goal of seliNDRIx is to construct the selection index from the
animal breeding data using the mixed and random sire models.

## Installation

You can install the development version of seliNDRIx like so:

\`\`\` r \# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
install.packages(“seliNDRIx”) require(seliNDRIx) \## Example This is a
basic example which shows you how to solve a common problem:

library(seliNDRIx)

# Example using mixed function

# Read the data

data <- data("data", package = "seliNDRIx") \# Define your parameters traits \<-
c(“tmy”, “py”, “fatyield”) fixed \<- c(“farm”, “soc”, “poc”) random \<-
c(“sire”) economic_values \<- c(1, 0.85, 0.65)

\# Run the analysis results \<- mixed_si( data = data, traits = traits,
fixed = fixed, random = random, economic_values = economic_values )
results \# To calculate the overall selection index for each animal SI
\<- c(results\$SelectionIndex) \# Selection index estimates (weights)
for traits traits \<- c(“tmy”, “py”, “fatyield”) \# Define the trait
columns to use overall_index \<- function(data, SI, traits) { \# Ensure
the number of weights matches the number of trait columns if (length(SI)
!= length(traits)) { stop(“The number of weights must match the number
of trait columns.”) } \# Select only the defined trait columns and
calculate the index data %\>% rowwise() %\>% mutate(Index =
sum(c_across(all_of(traits)) \* SI)) %\>% ungroup() } \# Calculate the
selection index result3 \<- overall_index(data, SI, traits) \# Print the
result print(result3)

\# Select the top 20% of animals with the highest selection index values
top20 \<- result3 %\>% arrange(desc(Index)) %\>% \# Sort by Index in
descending order slice_head(prop = 0.2) \# Select the top 20%

\# Example using
random function \# Read the data data \<- data("data", package = "seliNDRIx") \# Run
the analysis results2 \<- random_si(data, traits = c(“tmy”, “py”,
“fatyield”), economic_values = c(1, 0.85, 0.65)) results2 \# To
calculate the overall selection index for each animal SI \<-
c(results\$SelectionIndex) \# Selection index estimates (weights) for
traits traits \<- c(“tmy”, “py”, “fatyield”) \# Define the trait columns
to use overall_index \<- function(data, SI, traits) { \# Ensure the
number of weights matches the number of trait columns if (length(SI) !=
length(traits)) { stop(“The number of weights must match the number of
trait columns.”) } \# Select only the defined trait columns and
calculate the index data %\>% rowwise() %\>% mutate(Index =
sum(c_across(all_of(traits)) \* SI)) %\>% ungroup() }

\# Calculate the selection index result3 \<- overall_index(data, SI,
traits)

\# Print the result print(result3)

\# Select the top 20% of animals with the highest selection index values
top20 \<- result3 %\>% arrange(desc(Index)) %\>% \# Sort by Index in
descending order slice_head(prop = 0.2) \# Select the top 20%

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
