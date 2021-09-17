using Plots

struct asd
	a
end

function f()
	a = 2
	b = a[1]
	display(length(a))
end

f()
