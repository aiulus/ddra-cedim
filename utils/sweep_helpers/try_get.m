function v = try_get(f, g)
    try, v = f();
    catch, v = g();
    end
end