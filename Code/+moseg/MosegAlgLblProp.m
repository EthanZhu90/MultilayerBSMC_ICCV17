classdef MosegAlgLblProp < moseg.MosegAlgAffBase
    properties 
        Name = 'LabelPropagation';
        Options = moseg.MosegAlgLblPropOptions;
    end
    
    methods
        function obj = MosegAlgLblProp(startframe, varargin)
            obj@moseg.MosegAlgAffBase(startframe, varargin{:});
        end        
    end
        
    methods(Access = protected)
        function stepProp(obj)
            n = length(obj.CurrState.memTrajIds);
            
            if (n == 0)
                return;
            end

            % il - indices labeled current, iu - indices unlabeled current
            trajlen = obj.CurrState.memEndFrame - obj.CurrState.memStartFrame;
            il = find(obj.CurrState.memLbls ~= 0 & (trajlen >= min(obj.Options.conflen, obj.CurrState.frameNo - obj.StartFrame + 1)));
            iu = find(obj.CurrState.memLbls == 0 | (trajlen < min(obj.Options.conflen, obj.CurrState.frameNo - obj.StartFrame + 1)));
            
            % Do we have labeled data ? 
            if (isempty(il))              
                return;
            end
                        
            obj.CurrState.lowAff = find(sum(obj.CurrState.W) < eps);
            
            m = size(obj.CurrState.memProb,2);
            
            if (obj.Options.trackprob)
                Zl = obj.CurrState.memProb;
                Zl(iu,:) = 0;
            else
                % Transform labels to 1 hot
                Zl = zeros(n,m);                
                Zl(sub2ind([n m], il', obj.CurrState.memLbls(il)')) = 1;
                
                if (obj.Options.initprob)
                    Zl(iu,:) = obj.CurrState.memProb(iu,:);
                end
                
                % Uncomment to Initialize using nearest neighbors
                %[~, ix] = min(obj.CurrState.Dm(il,iu));                
                %Zl(iu,:) = obj.CurrState.memLbls(il(ix));                
            end
            
            if (obj.Options.lblpropglc)
                Znew = obj.lblpropglc(obj.CurrState.W, Zl);
            else
                if (obj.Options.fixlbl)
                    if (obj.Options.initprob)
                        Znew = obj.lblpropharmonic(obj.CurrState.W, Zl);
                        Znew(il,:) = 0;
                        Znew(il,1:end-1) = Zl(il,1:end);
                    else
                        Znew = obj.lblpropharmonic2(obj.CurrState.W, Zl, il, iu);
                    end
                else
                    Znew = obj.lblpropharmonic(obj.CurrState.W, Zl);
                end
            end

            % Make a probaility distribution
            % Avoid divide by zero problem
            ss = sum(Znew,2);
            ss(ss < eps) = 1; 
            Znew = bsxfun(@rdivide, Znew, ss);
            
            % Compute normalization
            if (obj.Options.normalize)
                if (size(Znew,2) > size(Zl,2))
                    eta = obj.Options.eta;
                    Zl = [(1-eta) * Zl eta * ones(size(Zl,1),1)];
                end
                pk = sum(Zl(il,:)) ./ length(il);
                mk = sum(Znew) ./ n; 
                mk(mk < eps & pk<eps) = 1; % No intial labels for this class
                mk(mk < eps) = 1; % No new labels for this class (ok let it dissapear.
                Znew = bsxfun(@rdivide, bsxfun(@times, Znew,pk), mk);
            end
            
            % Make decision
            [~,x] = max(Znew,[], 2);
            
            if (any(x == (m+1)))
                fprintf('****Introduced a new cluster******\n');
            else
                fprintf('\tNo new clusters introduced\n');
                Znew = Znew(:, 1:m);
            end                
            
            obj.CurrState.memLbls = x';
            obj.CurrState.memProb = Znew;             
        end
      
                      
        function Z = lblpropglc(obj, W, Zl)
                % Compute embedding using out-of-sampling.
                Z = moseg.MosegAlgLblProp.lblprop(W, Zl,obj.Options.alpha2,obj.Options.n_iter);            
        end
        
        function Z = lblpropharmonic(obj, W, Zl)
            n = size(Zl,1);
            eta = obj.Options.eta;
            eta2 = obj.Options.eta2; % Corresponds to eta in our paper
            Zlbar = [(1-eta)*Zl eta*ones(n,1)];
            %I = speye(n);
            ss = full(sum(W,2));
            ss(obj.CurrState.lowAff) = 1;
            Dinv = spdiags(1./ss, 0, n,n);
            P = Dinv * W;
            %Zu = (I - (1-eta)* P) \ Zlbar;
            %Zu = moseg.MosegUtils.solveiter(-(1-eta) .* P, I, Zlbar);
            Z = moseg.MosegUtils.solveiter2(-(1-eta) .* (1-eta2).* P, Zlbar);
        end
        
        function Z = lblpropharmonic2(obj, W, Zl, il, iu)
            n = size(Zl,1);
            eta = obj.Options.eta;
            %I = speye(n);
            ss = full(sum(W,2));
            ss(obj.CurrState.lowAff) = 1;
            Dinv = spdiags(1./ss, 0, n,n);
            P = Dinv * W;
            Puu = P(iu,iu);
            Pul = P(iu,il);
            b = [(1-eta) .* Pul * Zl(il,:) eta * ones(length(iu),1)];
            Zu = moseg.MosegUtils.solveiter2(-(1-eta) .* Puu, b);
            Z = [Zl zeros(n,1)];
            Z(iu,:) = Zu;
        end        
    end
    
    methods (Static)
        function Y = lblprop(W, Y, alpha, niter)
        % LBLPROP Label propagation on a Similarity graph.
 
            if (size(Y,1) == 1); Y =Y'; end;
            assert(all(all(W' == W)));

            n  = size(W,1);
            % Set W(i,i) to 0
            W(eye(n) == 1) = 0;

            % Compute Laplacian L
            Dm=diag(sum(W,2));
            D_p5=diag(diag(Dm).^(-.5));
            L=D_p5*W*D_p5;

            % Initialize Y0
            Y0= Y;
            % Iterate
            for j=1:1*niter
                Y=(1-alpha)*L*Y + alpha*Y0;              
            end
        end        
    end

end