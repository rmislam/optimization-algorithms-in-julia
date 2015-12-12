# First-order steepest descent potential reduction

include("steepest_descent.jl")

m = 8
n = 10
𝜀 = 0.5
max_iter = 20
converged_vec = zeros(max_iter)
reached_epsilon_vec = zeros(max_iter)
steps_vec = zeros(max_iter)
steps_limit_vec = zeros(max_iter)
cond_vec = zeros(max_iter)

for i = 1:max_iter
    println(i)
    converged_vec[i], reached_epsilon_vec[i], steps_vec[i], steps_limit_vec[i], cond_vec[i] = steepest_descent(m,n,𝜀)
end

writedlm("steepest_descent_results.txt", [converged_vec reached_epsilon_vec steps_vec steps_limit_vec cond_vec])


