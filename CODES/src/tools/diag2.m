function spec=diag2(A)
    for i=1:size(A,3)
       spec(:,i)=diag(A(:,:,i)); 
    end    
end