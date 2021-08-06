using MIPTM, Plots
function first_state(L, N, cap) #cap = max amount of bosons on one site...
    fullStates = Int(floor(N/cap))
    bosonsLeftOver = N % cap
    return [cap * ones(fullStates); bosonsLeftOver; zeros(L - fullStates - 1)]
end
function f()
	display(first_state(8, 8, 7))
end

function displayVec(v)
	i = q_find_index(v)
	display(i)
	show(v)
end

f()
