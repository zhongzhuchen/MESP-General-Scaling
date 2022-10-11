for s=2:62
    if fix63(s).nearzero_intgap == 0
        if ~isempty(intersect(fix63(s).fixto0list, fix63(s).fixto1list))
            fix63(s).fixto0list
            fix63(s).fixto1list
            s
            error("1");
        end
    end
end