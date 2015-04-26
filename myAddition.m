function myRes = myAddition( myOpr, var1, var2, var3 )

if (nargin==2)
	var3 = 100;
end

switch lower(myOpr)
	case {'add'}
		myRes = var1 + var2 + var3;
	case {'minus'}
		myRes = var1 - var2 - var3;
	case {'multi'}
		myRes = var1 * var2 * var3;
end
end

