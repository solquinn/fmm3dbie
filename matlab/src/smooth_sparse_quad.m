function Asmth = smooth_sparse_quad(kern,targs,S,row_ptr,col_ind,norderup,lbat)
    if nargin < 7
        lbat = 1e3;
    end
    if isa(kern,'kernel')
        kern = kern.eval;
    end

    nover = S.norders(1) + norderup; 
    S_over = oversample(S,nover);

    ixyzs = S_over.ixyzs;
    npols = ixyzs(2:end)-ixyzs(1:end-1);
    npts_col = npols(col_ind);

    [nt,~] = size(row_ptr);
    ntarg = nt-1;
    nrep = zeros(ntarg,1);
    istarts = row_ptr(1:end-1);
    iends = row_ptr(2:end)-1;
    icol_ind = zeros(sum(npts_col),1);
    isrcinds = cell(S_over.npatches,1);
    for i=1:S_over.npatches
        isrcinds{i} = ixyzs(i):ixyzs(i+1)-1;
    end

    istart = 1;
    for i=1:ntarg
        nrep(i) = sum(npts_col(istarts(i):iends(i))); 
        iinds = horzcat(isrcinds{col_ind(istarts(i):iends(i))});
        nelem = length(iinds);
        icol_ind(istart:istart+nelem-1) = iinds;
        istart = istart+nelem;
    end
    irow_ind = repelem((1:ntarg)',nrep);    
    
   nnz = length(irow_ind);
   vals = zeros(1,nnz);
    
    nbat = ceil(nnz/lbat);
    
    for k = 1:nbat
        ks = (lbat*(k-1)+1):min(lbat*k,nnz);
        rsrc = S_over.r(:,icol_ind(ks));
        rtarg = targs(:,irow_ind(ks));
        rnear = rtarg - rsrc;
        vals(ks) = kern(struct('r',[0;0;0]),struct('r',rnear)).*S_over.wts(icol_ind(ks));
    end
  
   Asmth = sparse(irow_ind,icol_ind, vals, ntarg, S_over.npts);

   rnodes = koorn.rv_nodes(S.norders(1));
    rnodes_over = koorn.rv_nodes(nover);

    vmat = koorn.coefs2vals(S.norders(1),rnodes_over);
    umat = koorn.vals2coefs(S.norders(1),rnodes);  
    Ainterp = vmat*umat;
        Ainterpc = cell(S.npatches,1);
    for i = 1:S.npatches
        Ainterpc{i} = sparse(Ainterp);
    end
    val2over = blkdiag(Ainterpc{:});
    Asmth = Asmth*val2over;
end