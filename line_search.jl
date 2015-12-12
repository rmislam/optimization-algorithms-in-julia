function line_search(x,d,A)
    𝛼_prev = 0
    𝛼_cur = 0.001
    𝛼_next = 𝛼_cur
    
    toler = 0.000001

    𝜀 = 0.2
    𝜂 = 2  # = 10
    𝜙𝛼 = (0.5 * (A*(x + 𝛼_next*d))' * (A*(x + 𝛼_next*d)))[1]  
    𝜙0 = (0.5 * (A*x)' * (A*x))[1]  
    d𝜙d0 = (d'*A'*A*x)[1] 
    original_satisfy = bool(𝜙𝛼 <= (𝜙0 + 𝜀*d𝜙d0*𝛼_next + toler))

    if bool(sum(x + 𝛼_cur*d .< toler))
        return 𝛼_prev
    end

    # These must be mutually exclusive
    use_quad = true
    use_cubic = false
    use_armijo = false

    armijo_pass = false

    steps = 0
    steps_limit = 10

    # Quadratic fit (method of false position)
    while !armijo_pass & (steps < steps_limit) & !use_armijo
        
        if use_quad
            dfd𝛼_cur = (d'*A'*A*(x + 𝛼_cur*d))[1]
            dfd𝛼_prev = (d'*A'*A*(x + 𝛼_prev*d))[1]
            𝛼_next = 𝛼_cur - dfd𝛼_cur*(𝛼_prev - 𝛼_cur)/(dfd𝛼_prev - dfd𝛼_cur)
        end

        if use_cubic 
            f𝛼_cur = (0.5 * (A*(x + 𝛼_cur*d))' * (A*(x + 𝛼_cur*d)))[1]  
            f𝛼_prev = (0.5 * (A*(x + 𝛼_prev*d))' * (A*(x + 𝛼_prev*d)))[1]  
            dfd𝛼_cur = (d'*A'*A*(x + 𝛼_cur*d))[1]
            dfd𝛼_prev = (d'*A'*A*(x + 𝛼_prev*d))[1]

            u1 = dfd𝛼_prev + dfd𝛼_cur - 3*(f𝛼_prev - f𝛼_cur)/(𝛼_prev - 𝛼_cur)
            u2_arg = u1^2 - dfd𝛼_prev*dfd𝛼_cur
            println("cubic")
            if u2_arg >= 0
                u2 = sqrt(u2_arg)
            else
                println("negative")
                break
            end
            𝛼_next = 𝛼_cur - (𝛼_cur - 𝛼_prev) * (dfd𝛼_cur + u2 - u1) / (dfd𝛼_cur - dfd𝛼_prev + 2*u2)
        end

        𝜀 = 0.2
        𝜂 = 2  # = 10
        𝜙𝛼 = (0.5 * (A*(x + 𝛼_next*d))' * (A*(x + 𝛼_next*d)))[1]  
        𝜙0 = (0.5 * (A*x)' * (A*x))[1]  
        d𝜙d0 = (d'*A'*A*x)[1] 
        𝜙𝜂𝛼 = (0.5 * (A*(x + 𝜂*𝛼_next*d))' * (A*(x + 𝜂*𝛼_next*d)))[1]  

        armijo_pass = (𝜙𝛼 <= (𝜙0 + 𝜀*d𝜙d0*𝛼_next + toler)) & (𝜙𝜂𝛼 > (𝜙0 + 𝜀*d𝜙d0*𝜂*𝛼_next - toler))

        if bool(sum(x + 𝛼_next*d .< toler))
            break
        end
        
        𝛼_prev = 𝛼_cur
        𝛼_cur = 𝛼_next
        steps += 1
    end
    
    while use_armijo
        𝜙𝛼 = (0.5 * (A*(x + 𝛼_next*d))' * (A*(x + 𝛼_next*d)))[1]  
        𝜙0 = (0.5 * (A*x)' * (A*x))[1]  
        d𝜙d0 = (d'*A'*A*x)[1] 

        if (𝜙𝛼 <= (𝜙0 + 𝜀*d𝜙d0*𝛼_next + toler)) & original_satisfy & !bool(sum(x + 𝛼_next*d .< toler))
            𝛼_cur = 𝛼_next
            𝛼_next *= 𝜂
        else
            break
        end
        if !original_satisfy  # no need to enforce nonnegativity constraint here since step size is shrinking
            if 𝜙𝛼 <= (𝜙0 + 𝜀*d𝜙d0*𝛼_next + toler)
                return 𝛼_next
            end
            𝛼_next /= 𝜂
        end
    end

    return 𝛼_cur
end

