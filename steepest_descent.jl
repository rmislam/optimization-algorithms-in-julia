# Homogeneous linear feasibility problem in Karmarkar form
#
#  A * x = 0
#  e' * x = 1
#  x >= 0
#
include("line_search.jl")
include("line_search_phi.jl")

function steepest_descent(m,n,𝜀)

    function update_conj(x,A,n,e)
        X = diagm(vec(x))
        𝛾 = 0.5 * maximum(eig(A'*A)[1])
        f = (0.5 * (A*x)'*(A*x))[1]
        ∇f = A'*A*x
        𝜌 = n + sqrt(n)
        ∇𝜙 = (𝜌/f)*∇f - inv(X) * e

        𝜆 = ((e' * X^2 * ∇𝜙 * f)/(𝜌 * x' * x))[1]
        p = X*((𝜌/f)*(∇f - e * 𝜆)) - e
        𝛽 = 1/(2 + 𝜌*𝛾/f)
        d_cur = -(𝛽/norm(p))*X*p
        g_cur = -d_cur

        angle = acos(dot(∇𝜙,d_cur)/(norm(∇𝜙)*norm(d_cur)))
        temp = readdlm("angle_measurements.txt")
        writedlm("angle_measurements.txt", [temp; angle])

        for k = 1:n
            𝛼 = line_search(x,d_cur,A)
            x = x + 𝛼*d_cur
            if k == n
                return x
            end

            X = diagm(vec(x))
            𝛾 = 0.5 * maximum(eig(A'*A)[1])
            f = (0.5 * (A*x)'*(A*x))[1]
            ∇f = A'*A*x
            𝜌 = n + sqrt(n)
            ∇𝜙 = (𝜌/f)*∇f - inv(X) * e

            𝜆 = ((e' * X^2 * ∇𝜙 * f)/(𝜌 * x' * x))[1]
            p = X*((𝜌/f)*(∇f - e * 𝜆)) - e
            𝛽 = 1/(2 + 𝜌*𝛾/f)
            d_add = -(𝛽/norm(p))*X*p
            g_next = -d_add

            angle = acos(dot(∇𝜙,d_add)/(norm(∇𝜙)*norm(d_add)))
            temp = readdlm("angle_measurements.txt")
            writedlm("angle_measurements.txt", [temp; angle])

            beta = ((g_next - g_cur)'*g_next/(g_cur'*g_cur))[1]

            d_cur = d_add + beta*d_cur
            g_cur = g_next
        end
    end

    function update(x,A,n,e)
        X = diagm(vec(x))
        𝛾 = 0.5 * maximum(eig(A'*A)[1])
        f = (0.5 * (A*x)'*(A*x))[1]
        ∇f = A'*A*x
        𝜌 = n + sqrt(n)
        ∇𝜙 = (𝜌/f)*∇f - inv(X) * e

        𝜆 = ((e' * X^2 * ∇𝜙 * f)/(𝜌 * x' * x))[1]
        p = X*((𝜌/f)*(∇f - e * 𝜆)) - e
        𝛽 = 1/(2 + 𝜌*𝛾/f)
        d = -(𝛽/norm(p))*X*p

        angle = acos(dot(∇𝜙,d)/(norm(∇𝜙)*norm(d)))
        temp = readdlm("angle_measurements.txt")
        writedlm("angle_measurements.txt", [temp; angle])

        # 𝛼 = 1  # if line search is not used
        𝛼 = line_search(x,d,A)

        return x + 𝛼*d
    end

    A = rand(-100:1:100,m,n)
    while (rank(A) < m) | !all(eigvals(A'*A) .>= 0) | isposdef(A'*A)
        A = rand(-100:1:100,m,n)
    end

    𝜌 = n + sqrt(n)
    e = ones(n)
    x0 = (1/n) * e

    x_cur = x0
    fx0 = (0.5 * (A*x0)' * (A*x0))[1]
    fx_cur = fx0
    fx_prev = fx0
    steps = 0

    𝛾 = 0.5 * maximum(eig(A'*A)[1])
    steps_limit = 4 * (n + sqrt(n)) * (1/𝜀) * log(1/𝜀) * maximum([1, (2*(n + sqrt(n))*𝛾)/(fx0[1]) ]) 

    converged = 0
    toler = 0.000001
    hard_limit = 1000

    while ((fx_cur/fx0)[1] > 𝜀 + toler) & (steps < hard_limit) 
        x_prev = x_cur
        x_cur = update_conj(x_prev,A,n,e)
        fx_prev = (0.5 * (A*x_prev)' * (A*x_prev))[1]
        fx_cur = (0.5 * (A*x_cur)' * (A*x_cur))[1]
        𝜙x_prev = 𝜌*log(fx_prev) - sum(log(x_prev))
        𝜙x_cur = 𝜌*log(fx_cur) - sum(log(x_cur))
        steps += 1
        if 𝜙x_prev - 𝜙x_cur < toler
        #if fx_prev - fx_cur < toler
            converged = 1
            break
        end
    end 

    reached_epsilon = 1
    if (fx_cur/fx0)[1] > 𝜀 + toler 
        reached_epsilon = 0
    end

    return converged, reached_epsilon, steps, steps_limit, cond(A)
end


