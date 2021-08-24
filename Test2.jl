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
	y = 3
	for i in 1:3
		local x = 2
	end
	display(x)
end

f()
