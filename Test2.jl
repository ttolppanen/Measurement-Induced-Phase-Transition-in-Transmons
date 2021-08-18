using ParametersModule, MIPTM

struct unmutableStruct
	a
	b
	function unmutableStruct(x)
		new(x, 2*x)
	end
end
mutable struct mutableStruct
	c
	d
end

function f()
	a = [1 2 3; 3 3 3]
	b = a
	b .+= [1 1 1; 4 4 4]
	display(b)
	display(a)
end

f()
