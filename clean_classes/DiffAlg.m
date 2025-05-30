classdef DiffAlg < handle
    % DiffAlg implements a truncated power series (differential algebra)
    % in n variables (nDv), following the framework introduced by Berz.
    %
    % The algebra is defined by a power series expansion truncated at a given
    % total degree (Order). The coefficients are stored in a vector (Series)
    % whose length is given by:
    %
    %       dim = nchoosek(nVars + Order, nVars)
    %
    % The multi-index bookkeeping is handled via the getMultiIndices method,
    % which returns a matrix where each row is an exponent vector corresponding
    % to a term in the series.
    %
    % CONSTRUCTOR USAGE:
    %
    % 1. Explicit series:
    %    DA = DiffAlg(series, Order, nVars);
    %
    % 2. Construct a DA object representing a specific variable:
    %    DA = DiffAlg('var', varIndex, Order, nVars);
    %
    % 3. Construct a DA object representing the constant 1:
    %    DA = DiffAlg('one', Order, nVars);
    %
    % Alternatively, use the static methods DiffAlg.var(varIndex, Order, nVars)
    % and DiffAlg.one(nVars, Order).
    %
    % Additional operations include addition, subtraction, multiplication,
    % in-place multiplication, differentiation, exponentiation, and now division.
    
    properties
        Series  % Coefficient vector for the truncated power series.
        Order   % Truncation order (maximum total degree).
        nVars   % Number of variables in the power series.
    end
    
    properties (Access = private)
        MultiIndicesCache  % Cache for multi-index bookkeeping.
    end
    
    methods
        function result = evaluate(obj, initial_values, max_eval_order)
            % Evaluates the DA object (power series) for given initial values,
            % optionally truncating the evaluation at a specified maximum order.
            % INPUTS:
            %   obj: A DiffAlg object.
            %   initial_values: A row vector [value_var1, value_var2, ...]
            %   max_eval_order: (Optional) Maximum total order to include
            %                   in the evaluation. Defaults to obj.Order.
            % OUTPUTS:
            %   result: The scalar numerical value of the series evaluation.

            arguments
                obj DiffAlg
                initial_values (1,:) double
                max_eval_order (1,1) double = obj.Order % Default to full order
            end

            if length(initial_values) ~= obj.nVars
                error('DiffAlg:evaluate:InputSizeMismatch', ...
                      'Number of initial values (%d) must match number of DA variables (%d).', ...
                      length(initial_values), obj.nVars);
            end
             if max_eval_order < 0
                 warning('DiffAlg:evaluate:NegativeOrder', 'max_eval_order cannot be negative. Using 0.');
                 max_eval_order = 0;
             elseif max_eval_order > obj.Order
                  % It's okay to request a higher order, just evaluate up to obj.Order
                  max_eval_order = obj.Order;
             end
            if isempty(obj.Series)
                warning('DiffAlg:evaluate:EmptySeries', 'Attempting to evaluate an empty DA object. Returning 0.');
                result = 0.0;
                return;
            end

            coeffs = obj.Series;
            MI = obj.getMultiIndices(); % Multi-indices matrix
            result = 0.0; % Use double precision

            if size(MI, 1) ~= length(coeffs)
                error('DiffAlg:evaluate:InternalInconsistency', ...
                      'Internal mismatch between number of coefficients (%d) and multi-indices (%d).', ...
                      length(coeffs), size(MI, 1));
            end

            % Calculate total degree for each term (only need to do this once)
            term_orders = sum(MI, 2);

            for k = 1:length(coeffs)
                % --- Check evaluation order ---
                if term_orders(k) > max_eval_order
                    continue; % Skip this term if its order is too high
                end

                term_coeff = coeffs(k);
                if term_coeff == 0
                    continue;
                end

                monomial_value = 1.0;
                exponents = MI(k, :);

                for j = 1:obj.nVars
                    exp_j = exponents(j);
                    if exp_j == 0, continue; end % x^0 = 1
                    if exp_j == 1
                        monomial_value = monomial_value * initial_values(j);
                    else % exp_j > 1
                        if initial_values(j) == 0
                            monomial_value = 0.0;
                            break;
                        end
                        monomial_value = monomial_value * (initial_values(j) ^ exp_j);
                    end
                    if ~isfinite(monomial_value), break; end
                end

                 term_value = term_coeff * monomial_value;
                 if ~isfinite(term_value)
                     warning('DiffAlg:evaluate:NonFiniteTerm', ...
                             'Non-finite term value encountered for term k=%d, order=%d at eval_order=%d. Skipping term.', k, term_orders(k), max_eval_order);
                     continue; % Skip non-finite term
                 end
                result = result + term_value;
                 if ~isfinite(result)
                     warning('DiffAlg:evaluate:NonFiniteSum', ...
                             'Sum became non-finite after adding term k=%d. Returning current result.', k);
                     return; % Return early if sum becomes non-finite
                 end
            end
        end

        %% Constructor with multiple options:
        %   1. DiffAlg(series, Order, nVars)
        %   2. DiffAlg('var', varIndex, Order, nVars)
        %   3. DiffAlg('one', Order, nVars)
        function obj = DiffAlg(varargin)
            if nargin == 0
                obj.Series = [];
                obj.Order = 0;
                obj.nVars = 0;
                obj.MultiIndicesCache = [];
            elseif ischar(varargin{1})
                switch lower(varargin{1})
                    case 'var'
                        if nargin < 4
                            error('Usage: DiffAlg(''var'', varIndex, Order, nVars)');
                        end
                        varIndex = varargin{2};
                        order = varargin{3};
                        nVars = varargin{4};
                        expectedDim = nchoosek(nVars + order, nVars);
                        series = zeros(expectedDim, 1);
                        tempDA = DiffAlg(series, order, nVars);
                        MI = tempDA.getMultiIndices();
                        target = zeros(1, nVars);
                        target(varIndex) = 1;
                        idx = find(all(MI == target, 2), 1);
                        if isempty(idx)
                            error('Multi-index for variable %d not found.', varIndex);
                        end
                        series(idx) = 1;
                        obj.Series = series;
                        obj.Order = order;
                        obj.nVars = nVars;
                        obj.MultiIndicesCache = [];
                    case 'one'
                        if nargin < 3
                            error('Usage: DiffAlg(''one'', Order, nVars)');
                        end
                        order = varargin{2};
                        nVars = varargin{3};
                        expectedDim = nchoosek(nVars + order, nVars);
                        series = zeros(expectedDim, 1);
                        series(1) = 1;
                        obj.Series = series;
                        obj.Order = order;
                        obj.nVars = nVars;
                        obj.MultiIndicesCache = [];
                    otherwise
                        error('Unknown constructor option: %s', varargin{1});
                end
            else
                series = varargin{1};
                order = varargin{2};
                nVars = varargin{3};
                expectedDim = nchoosek(nVars + order, nVars);
                if isempty(series)
                    series = zeros(expectedDim, 1);
                elseif numel(series) ~= expectedDim
                    error('The series length must be (nVars+Order choose nVars) = %d.', expectedDim);
                end
                obj.Series = series;
                obj.Order = order;
                obj.nVars = nVars;
                obj.MultiIndicesCache = [];
            end
        end
        
        %% getMultiIndices returns the multi-index matrix.
        function indices = getMultiIndices(obj)
            if isempty(obj.MultiIndicesCache)
                obj.MultiIndicesCache = DiffAlg.generateMultiIndices(obj.Order, obj.nVars);
            end
            indices = obj.MultiIndicesCache;
        end
        
        %% Overload plus operator.
        function result = plus(a, b)
            if a.Order ~= b.Order || a.nVars ~= b.nVars
                error('Both series must have the same Order and nVars.');
            end
            newSeries = a.Series + b.Series;
            result = DiffAlg(newSeries, a.Order, a.nVars);
        end
        
        %% Overload minus operator.
        function result = minus(a, b)
            if a.Order ~= b.Order || a.nVars ~= b.nVars
                error('Both series must have the same Order and nVars.');
            end
            newSeries = a.Series - b.Series;
            result = DiffAlg(newSeries, a.Order, a.nVars);
        end
        
        %% Overload mtimes operator for multiplication.
        function result = mtimes(a, b)
            if isa(a, 'DiffAlg') && isa(b, 'DiffAlg')
                if a.Order ~= b.Order || a.nVars ~= b.nVars
                    error('DA objects must have the same Order and nVars for multiplication.');
                end
                order = a.Order;
                nVars = a.nVars;
                MI = a.getMultiIndices();
                N = size(MI, 1);
                newSeries = zeros(N, 1);
                for i = 1:N
                    for j = 1:N
                        gamma = MI(i,:) + MI(j,:);
                        if sum(gamma) <= order
                            k = find(all(MI == gamma, 2), 1);
                            if ~isempty(k)
                                newSeries(k) = newSeries(k) + a.Series(i) * b.Series(j);
                            end
                        end
                    end
                end
                result = DiffAlg(newSeries, order, nVars);
            elseif isa(a, 'DiffAlg') && isa(b, 'double')
                result = DiffAlg(a.Series * b, a.Order, a.nVars);
            elseif isa(a, 'double') && isa(b, 'DiffAlg')
                result = DiffAlg(a * b.Series, b.Order, b.nVars);
            else
                error('Multiplication not defined for these operand types.');
            end
        end
        
        %% In-place multiplication.
        function multiplyAssign(obj, operand)
            if isa(operand, 'DiffAlg')
                if obj.Order ~= operand.Order || obj.nVars ~= operand.nVars
                    error('DA objects must have the same Order and nVars for multiplication.');
                end
                order = obj.Order;
                nVars = obj.nVars;
                MI = obj.getMultiIndices();
                N = size(MI, 1);
                newSeries = zeros(N, 1);
                for i = 1:N
                    for j = 1:N
                        gamma = MI(i,:) + MI(j,:);
                        if sum(gamma) <= order
                            k = find(all(MI == gamma, 2), 1);
                            if ~isempty(k)
                                newSeries(k) = newSeries(k) + obj.Series(i) * operand.Series(j);
                            end
                        end
                    end
                end
                obj.Series = newSeries;
            elseif isa(operand, 'double')
                obj.Series = obj.Series * operand;
            else
                error('Unsupported operand type for in-place multiplication.');
            end
        end
        
        %% differentiate computes the derivative with respect to variable varIndex.
        % The result is a DA object of the same Order, with zeros filling in the
        % terms that are truncated.
        function deriv = differentiate(obj, varIndex)
            if varIndex < 1 || varIndex > obj.nVars
                error('varIndex must be between 1 and nVars.');
            end
            newOrder = obj.Order;
            newDim = nchoosek(obj.nVars + newOrder, obj.nVars);
            newSeries = zeros(newDim, 1);
            MI = obj.getMultiIndices();
            % For each term in the original series, if the exponent for varIndex is > 0,
            % add its contribution to the term corresponding to (alpha - e_i)
            for i = 1:length(obj.Series)
                currentCoeff = obj.Series(i);
                alpha = MI(i,:);
                if alpha(varIndex) > 0
                    new_alpha = alpha;
                    new_alpha(varIndex) = new_alpha(varIndex) - 1;
                    k = find(all(MI == new_alpha, 2), 1);
                    if ~isempty(k)
                        newSeries(k) = newSeries(k) + currentCoeff * alpha(varIndex);
                    end
                end
            end
            deriv = DiffAlg(newSeries, newOrder, obj.nVars);
        end
        
        %% exp computes the exponential of the DA object.
        function result = exp(obj)
            kMax = floor(obj.Order/2) + 1;
            result = DiffAlg.one(obj.nVars, obj.Order);
            term = DiffAlg.one(obj.nVars, obj.Order);
            for k = 1:kMax
                term = term * obj;
                result = result + term * (1/factorial(k));
            end
        end
        
        %% inverse computes the multiplicative inverse of the DA object.
        % That is, for F such that F*G = 1, this method returns G.
        function invObj = inverse(obj)
            % The constant term (degree 0) is the first element in the Series.
            a0 = obj.Series(1);
            if abs(a0) < eps
                error('Multiplicative inverse does not exist: zero constant term.');
            end
            N = numel(obj.Series);
            newSeries = zeros(N,1);
            newSeries(1) = 1 / a0;
            MI = obj.getMultiIndices();  % Assume sorted in increasing total degree
            % Loop over all terms i > 1
            for i = 2:N
                beta = MI(i, :);
                sum_term = 0;
                % Sum over all terms j (excluding constant term) that can contribute
                % to beta. (We assume that MI is sorted so that terms with lower degree come first.)
                for j = 2:N
                    alpha = MI(j, :);
                    if all(beta >= alpha)
                        diff = beta - alpha;
                        k = find(all(MI == diff, 2), 1);
                        if ~isempty(k)
                            sum_term = sum_term + obj.Series(j) * newSeries(k);
                        end
                    end
                end
                newSeries(i) = - (1 / a0) * sum_term;
            end
            invObj = DiffAlg(newSeries, obj.Order, obj.nVars);
        end
        
        %% Overload division operator (mrdivide) so that A/B = A * inverse(B).
        function result = mrdivide(a, b)
            if isa(b, 'DiffAlg')
                result = a * b.inverse();
            elseif isa(a, 'DiffAlg') && isa(b, 'double')
                result = DiffAlg(a.Series / b, a.Order, a.nVars);
            else
                error('Division is only defined for DA objects or DA divided by scalar.');
            end
        end

        function sqrtObj = dsqrt(obj)
            % dsqrt computes the truncated square-root expansion of a DA object 'obj'.
            % The constant term must be > 0 to avoid singularities.
            
            a0 = obj.Series(1);
            if a0 <= 0.0
                error('Square root does not exist: the constant term is non-positive.');
            end
            
            N = numel(obj.Series);
            newSeries = zeros(N,1);
            
            % constant term is sqrt(a0)
            s0 = sqrt(a0);
            newSeries(1) = s0;
            
            % We'll need the forward approach. Suppose F(x)=a0 + ...,
            % and G(x)=sqrt(F(x)) = s0 + ...
            % Then for each higher-order term, we do:
            %   g_beta = (1/(2 s0)) [ f_beta - sum(...) ],
            % similar to the multiplicative inverse expansions.
            
            MI = obj.getMultiIndices();
            for i = 2:N
                beta = MI(i,:);
                % sum_term = sum_{alpha>0, alpha<beta} f_alpha*g_{beta-alpha}
                % (like a convolution). We'll gather that from the existing newSeries.
                sum_contrib = 0.0;
                for j = 2:N
                    alpha = MI(j,:);
                    if all(beta >= alpha) && any(alpha>0)
                        diffIdx = beta - alpha;
                        k = find(all(MI == diffIdx, 2), 1);
                        if ~isempty(k)
                            sum_contrib = sum_contrib + obj.Series(j)*newSeries(k);
                        end
                    end
                end
                % then g_beta = (1/(2 s0)) * [ f_beta - sum_contrib ]
                fb = obj.Series(i);
                newSeries(i) = (1.0/(2.0*s0))*(fb - sum_contrib);
            end
            
            sqrtObj = DiffAlg(newSeries, obj.Order, obj.nVars);
        end

        function displaySeries(obj, precision, varNames, threshold)
            % displaySeries Displays the DA series coefficients grouped by order.
            %
            % INPUTS:
            %   obj: The DiffAlg object.
            %   precision (optional): Format specifier string for coefficients
            %                      (e.g., '%.4e', '%.6f', default '%.4e').
            %   varNames (optional): Cell array of strings for variable names
            %                      (e.g., {'r', 'p_r'}). Defaults to {'x1', 'x2', ...}.
            %   threshold (optional): Minimum absolute coefficient value to display.
            %                       Defaults to MATLAB's eps. Set to 0 to display all.

            % --- Input Handling ---
            if nargin < 2 || isempty(precision)
                precision = '%.4e'; % Default format
            end

            numVars = obj.nVars; % If obj is not a DiffAlg object, error occurs here.

            if nargin < 3 || isempty(varNames)
                varNames = arrayfun(@(i) sprintf('x%d', i), 1:numVars, 'UniformOutput', false);
            elseif ~iscell(varNames) || numel(varNames) ~= numVars 
                error('varNames must be a cell array with nVars (%d) elements.', numVars);
            end
             if nargin < 4 || isempty(threshold)
                threshold = eps; % Default threshold
            end

            % --- Get Data ---
            if isempty(obj.Series)
                fprintf('DA object is empty.\n');
                return;
            end
            MI = obj.getMultiIndices(); % If obj is not valid, might error here too.
            coeffs = obj.Series;
            Order_loc = obj.Order;
            nVars_loc = numVars;

            % --- Print Header ---
            fprintf('Displaying DA Series (Order=%d, nVars=%d)\n', Order_loc, nVars_loc);
            fprintf('%s\n', repmat('-', 1, 40)); % Separator line

            % --- Loop by Order ---
            % Calculate max length for index string alignment based on Order and nVars
            test_alpha = repmat(Order_loc, 1, nVars_loc); % Hypothetical worst-case index
            indexStr_example = sprintf('%d,', test_alpha);
            indexStr_example = sprintf('(%s)', indexStr_example(1:end-1));
            maxIndexStrLen = length(indexStr_example) + 2; % Add some padding

            maxCoeffStrLen = 0; % Find max length for coefficient alignment

            % Pre-calculate coefficient strings to determine padding
            coeffStrs = cell(numel(coeffs), 1);
             for i = 1:numel(coeffs)
                 if abs(coeffs(i)) >= threshold || threshold == 0
                      coeffStrs{i} = sprintf(precision, coeffs(i));
                      maxCoeffStrLen = max(maxCoeffStrLen, length(coeffStrs{i}));
                 else
                     coeffStrs{i} = ''; % Mark as below threshold
                 end
             end
             maxCoeffStrLen = max(maxCoeffStrLen, length('-Inf')); % Ensure space for potential Inf/NaN

            for k = 0:Order_loc
                fprintf('--- Order %d ---\n', k);
                indices_k = find(sum(MI, 2) == k); % Find indices for this order
                found_term_in_order = false;

                if isempty(indices_k) && k==0 && ~isempty(coeffs) % Handle case where only order 0 exists
                     indices_k = find(sum(MI,2)==0); % Should be index 1
                elseif isempty(indices_k)
                     % No terms for this order > 0, print message later if needed
                end


                num_terms_in_order = numel(indices_k);
                num_displayed_in_order = 0;

                for i = 1:num_terms_in_order
                    idx = indices_k(i);
                     % Ensure idx is valid for coeffs and coeffStrs
                     if idx > numel(coeffs) || idx < 1
                         warning('Invalid index encountered for order %d. Skipping.', k);
                         continue;
                     end

                    c = coeffs(idx);
                    coeffStr = coeffStrs{idx};

                    % Skip if below threshold (and threshold > 0)
                    if isempty(coeffStr) && threshold > 0
                        continue;
                    elseif isempty(coeffStr) && threshold == 0 % Means coeff was exactly 0.0
                         coeffStr = sprintf(precision, 0.0); % Display it if threshold is 0
                         maxCoeffStrLen = max(maxCoeffStrLen, length(coeffStr)); % Update max length if needed
                    end

                    num_displayed_in_order = num_displayed_in_order + 1;
                    found_term_in_order = true; % Mark that we are displaying something

                    alpha = MI(idx, :);

                    % Format Multi-index String: e.g., "(1,0,2)"
                    indexStr = sprintf('%d,', alpha);
                    indexStr = sprintf('(%s)', indexStr(1:end-1)); % Remove trailing comma

                    % Format Monomial String: e.g., "[x1 * x3^2]"
                    monoParts = {};
                    isConstant = true;
                    for j = 1:nVars_loc
                        if alpha(j) == 1
                            monoParts{end+1} = varNames{j};
                            isConstant = false;
                        elseif alpha(j) > 1
                            monoParts{end+1} = sprintf('%s^%d', varNames{j}, alpha(j));
                            isConstant = false;
                        end
                    end

                    if isConstant
                        monoStr = '[Constant]';
                    else
                        monoStr = ['[' strjoin(monoParts, ' * ') ']'];
                    end

                    % Print Line with padding for alignment
                    fprintf('  %-*s: %*s   %s\n', ...
                           maxIndexStrLen, indexStr, ... % Padded index
                           maxCoeffStrLen, coeffStr, ... % Padded coefficient
                           monoStr);                    % Monomial string
                end

                % Check if we displayed anything for this order
                if num_displayed_in_order == 0 && num_terms_in_order > 0 % Only print if terms existed but were filtered
                    if threshold > 0
                        fprintf('  (All terms zero or below threshold=%.2e)\n', threshold);
                    else % threshold == 0
                        fprintf('  (All terms zero)\n');
                    end
                elseif num_displayed_in_order == 0 && num_terms_in_order == 0 && k > 0
                     % If no terms were even defined for this order (k>0), print nothing
                     % fprintf('  (No terms defined for this order)\n'); % Optional, can be verbose
                end
            end

            % --- Print Footer ---
            fprintf('%s\n', repmat('-', 1, 40)); % Separator line
            if threshold > 0
                 fprintf('(Coefficients with abs value < %.2e not shown)\n', threshold);
            else
                 fprintf('(Showing all coefficients)\n');
            end
        end
    end
    
    methods (Static)
        %% one returns a DA object representing the constant 1.
        function obj = one(nVars, order)
            dim = nchoosek(nVars + order, nVars);
            series = zeros(dim, 1);
            series(1) = 1;
            obj = DiffAlg(series, order, nVars);
        end
        
        %% var returns a DA object representing the specified variable.
        % Usage: obj = DiffAlg.var(varIndex, Order, nVars);
        function obj = var(varIndex, order, nVars)
            if varIndex < 1 || varIndex > nVars
                error('varIndex must be between 1 and nVars.');
            end
            dim = nchoosek(nVars + order, nVars);
            series = zeros(dim, 1);
            temp = DiffAlg(series, order, nVars);
            MI = temp.getMultiIndices();
            target = zeros(1, nVars);
            target(varIndex) = 1;
            idx = find(all(MI == target, 2), 1);
            if isempty(idx)
                error('Multi-index for variable %d not found.', varIndex);
            end
            series(idx) = 1;
            obj = DiffAlg(series, order, nVars);
        end
    end
    
    methods (Static, Access = private)
        %% generateMultiIndices generates all multi-indices for given nVars and order.
        function indices = generateMultiIndices(order, nVars)
            indices = [];
            for total = 0:order
                comp = DiffAlg.compositions(total, nVars);
                indices = [indices; comp];  %#ok<AGROW>
            end
        end
        
        %% compositions recursively generates all compositions of n into k parts.
        function comps = compositions(n, k)
            if k == 1
                comps = n;
            else
                comps = [];
                for i = 0:n
                    subcomps = DiffAlg.compositions(n - i, k - 1);
                    comps = [comps; [repmat(i, size(subcomps, 1), 1), subcomps]]; %#ok<AGROW>
                end
            end
        end
    end
end
