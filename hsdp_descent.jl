# Homogeneous and self-dual linear feasibility problem
#
#

include("line_search_hsdp_x.jl")
include("line_search_hsdp_y.jl")

function hsdp_descent(m,n,𝜀)

    function update_x(x,y,A,B,n,e)
        X = diagm(vec(x))
        𝛾 = 0.5 * maximum(eig([A B]'*[A B])[1])
        f = (0.5 * ([A B]*[x; y])'*([A B]*[x; y]))[1]
        ∇fx = A'*A*x + A'*B*y
        𝜌 = n + sqrt(n)
        ∇𝜙x = (𝜌/f)*∇fx - inv(X) * e

        𝜆 = ((e' * X^2 * ∇𝜙x * f)/(𝜌 * x' * x))[1]
        p = X*((𝜌/f)*(∇fx - e * 𝜆)) - e
        𝛽 = 1/(2 + 𝜌*𝛾/f)
        dx = -(𝛽/norm(p))*X*p

        angle = acos(dot(∇𝜙x,dx)/(norm(∇𝜙x)*norm(dx)))
        temp = readdlm("angle_measurements_hsdp.txt")
        writedlm("angle_measurements_hsdp.txt", [temp; angle])

        𝛼 = 1  # if line search is not used
        # 𝛼 = line_search_hsdp_x(x,y,dx,A,B)

        return x + 𝛼*dx
    end

    function update_y(x,y,A,B)
        # unconstrained gradient descent
        𝛾 = 0.5 * maximum(eig([A B]'*[A B])[1])
        ∇fy = B'*A*x + B'*B*y
        dy = -(1/𝛾)*∇fy

        𝛼 = 1  # if line search is not used
        # 𝛼 = line_search_hsdp_y(x,y,dy,A,B)

        return y + 𝛼*dy
    end

    A0 = rand(-100:1:100,m,n)
    b = rand(-100:1:100,m)
    c = rand(-100:1:100,n)
    I = eye(n)
    A = [A0 zeros(m,n) -b zeros(m,1); zeros(n,n) -I c zeros(n,1); -c' zeros(1,n) 0 -1]
    B = [zeros(m,m); -A0'; b']
    
    while (rank(A) < m) | !all(eigvals(A'*A) .>= 0) | isposdef(A'*A)
        A0 = rand(-100:1:100,m,n)
        b = rand(-100:1:100,m)
        c = rand(-100:1:100,n)
        I = eye(n)
        A = [A0 zeros(m,n) -b zeros(m,1); zeros(n,n) -I c zeros(n,1); -c' zeros(1,n) 0 -1]
        B = [zeros(m,m); -A0'; b']
    end

    n = size(A,2)   # redefine n
    𝜌 = n + sqrt(n)
    e = ones(n)
    x0 = (1/n) * e
    y0 = (1/m) * ones(m) # y0 (initial y) is arbitrary. May want to choose a smarter initial point
    
    x_cur = x0
    y_cur = y0
    f0 = (0.5 * ([A B]*[x0; y0])' * ([A B]*[x0; y0]))[1]
    f_cur = f0
    f_prev = f0
    steps = 0

    𝛾 = 0.5 * maximum(eig([A B]'*[A B])[1])
    steps_limit = 4 * (n + sqrt(n)) * (1/𝜀) * log(1/𝜀) * maximum([1, (2*(n + sqrt(n))*𝛾)/(f0[1]) ]) 

    converged = 0
    toler = 0.000001
    hard_limit = 10000

    while ((f_cur/f0)[1] > 𝜀 + toler) & (steps < hard_limit) 
        x_prev = x_cur
        y_prev = y_cur
        x_cur = update_x(x_prev,y_prev,A,B,n,e)
        y_cur = update_y(x_prev,y_prev,A,B)
        f_prev = (0.5 * ([A B]*[x_prev; y_prev])' * ([A B]*[x_prev; y_prev]))[1]
        f_cur = (0.5 * ([A B]*[x_cur; y_cur])' * ([A B]*[x_cur; y_cur]))[1]
        𝜙x_prev = 𝜌*log(f_prev) - sum(log(x_prev))
        𝜙x_cur = 𝜌*log(f_cur) - sum(log(x_cur))
        steps += 1
        if 𝜙x_prev - 𝜙x_cur < toler
        #if f_prev - f_cur < toler
            converged = 1
            break
        end
    end 

    reached_epsilon = 1
    if (f_cur/f0)[1] > 𝜀 + toler 
        reached_epsilon = 0
    end

    return converged, reached_epsilon, steps, steps_limit, cond([A B])
end

