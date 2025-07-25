function downsampled = f_downsample(sig,bin)

if bin == 1
    downsampled = sig;
else   
    dim = size(sig);
    if size(dim,2) == 2
        dim(3) = 1;
    end
    downsampled = zeros(floor(dim(1)/bin),floor(dim(2)/bin),dim(3));
    
    for h = bin:bin:dim(1)
        for w = bin:bin:dim(2)
            downsampled(h/bin,w/bin,:) = mean(sig((h-bin+1):h,(w-bin+1):w,:),[1,2],'omitnan');
        end
    end

    check = zeros(floor(dim(1)/bin),floor(dim(2)/bin),dim(3));
    check_sig = ~isnan(mean(sig,3));
    for h = bin:bin:dim(1)
        for w = bin:bin:dim(2)
            check(h/bin,w/bin,:) = sum(check_sig((h-bin+1):h,(w-bin+1):w,:),[1,2]);
        end
    end
    check = check/bin.^2;
    check = double(check > 0.5);
    check(check==0) = NaN;
    downsampled = downsampled.*check;

end
