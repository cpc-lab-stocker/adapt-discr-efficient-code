function [ x1_straight, probCmpCw, bias, thresh ] = ECAdapt_2AFC_par( k_o, theta_o, kappa_o, k, theta, kappa, kappa_i, kappa_e, res, adaptor, stimulus )
% 2AFC: theta vs theta+dtheta
% stimiulus = test + adaptor; test must have res resolution

F_o = @(x)sum_n_vmpdf_180( x, k_o, theta_o, kappa_o );
ada_trans = integral(F_o, 0, adaptor) * 180; % adaptor in internal space w/o adaptation

x1 = circ180(adaptor-90:res:adaptor+90-res); % physical space, [0,180), adaptor in middle
x1_straight = circshift(x1, round((adaptor+90)/res)); % 0~180 in increasing order
x2 = NaN(1, length(x1));
parfor i = 1:length(x2)
    x2(i) = circ90(integral(F_o, 0, x1(i)) * 180 - ada_trans);  % internal space w/o adaptation, 0 = adaptor
end
prior0 = F_o(x1);
% dx1 = ones(1, length(x1))*res;
prior2 = 1/180;
f1 = prior0/prior2;

%% adapt
prior1 = NaN(1, length(x2));
x31 = NaN(1, length(x2));
F = @(x)sum_n_vmpdf_180( x, k, theta, kappa );
parfor i = 1:length(x2)
    cdf1 = integral(F,-90,x2(i)); 
    x31(i) = cdf1*180-90; % internal space w/ adaptation, 0 = adaptor
    prior1(i) = sum_n_vmpdf_180( x2(i), k, theta, kappa ); 
end
f21 = prior1/prior2;

p1_m1_theta = zeros(length(x1), length(x1)); % P( external sample m1 | stimulus theta ); m1 * theta
if kappa_e > 700
    p1_m1_theta = eye(length(x1));
else
    parfor j = 1:length(x1)
        p1_m1_theta(:,j) = circ_vmpdf(2*x1/180*pi,2*x1(j)/180*pi,kappa_e)'/90*pi;
    end
end
p1_m2_m1 = zeros(length(x1), length(x1)); % P( internal sample m2 | external sample m1 ); m2 * m1
if kappa_i > 700
    p1_m2_m1 = eye(length(x1));
else
    parfor j = 1:length(x1)
        % uniform von mises noise in sensory space, then transformed to physical space
        p1_m2_m1(:,j) = circ_vmpdf(2*x31/180*pi,2*x31(j)/180*pi,kappa_i)'/90*pi.*f21'.*f1';
    end
end
p1_m2_theta = p1_m2_m1 * p1_m1_theta;
p1_m2_theta = p1_m2_theta./(sum(p1_m2_theta,1)*res);

%% 2AFC
probCmpCw = NaN(length(stimulus), length(x1));
% probCmpCcw = NaN(length(stimulus), length(x1));
L = length(x1);
is_cw = [1/2, ones(1,L/2-1), 1/2, zeros(1,L/2-1)];
is_cw = circshift(flip(is_cw),1); % correct order for cconv
is_ccw = 1-is_cw;
parfor i = 1:length(stimulus)
    sti = stimulus(i);
    sti_ind = round(circ180(sti-adaptor+90)/res)+1;
    for j = 1:L %comparison stimulus
        p_cmp_cw = cconv(p1_m2_theta(:,j)', is_cw, L) * res;
        p_cmp_ccw = cconv(p1_m2_theta(:,j)', is_ccw, L) * res;
        
        prob_cmp_cw = p_cmp_cw * p1_m2_theta(:,sti_ind) * res;
        prob_cmp_ccw = p_cmp_ccw * p1_m2_theta(:,sti_ind) * res;
        
        probCmpCw(i,j) = prob_cmp_cw/(prob_cmp_cw+prob_cmp_ccw);
%         probCmpCcw(i,j) = prob_cmp_ccw/(prob_cmp_cw+prob_cmp_ccw);
    end
end

probCmpCw = circshift(probCmpCw, round((adaptor+90)/res), 2); % stimulus * x1_straight
% probCmpCcw = circshift(probCmpCcw, round((adaptor+90)/res), 2);

%% bias & threshold
if nargout > 2
    POE1 = NaN(1,length(stimulus));
    thresh = NaN(1,length(stimulus));

    [min_prob,min_prob_ind] = min(probCmpCw, [], 2);
    [max_prob,max_prob_ind] = max(probCmpCw, [], 2);

    for i = 1:length(stimulus)
        % put prob in increasing order
        if min_prob_ind(i) < max_prob_ind(i)
            prob_temp = probCmpCw(i, min_prob_ind(i):max_prob_ind(i));
            x1_temp = x1_straight(min_prob_ind(i):max_prob_ind(i));
        else
            prob_temp = probCmpCw(i, [min_prob_ind(i):size(probCmpCw,2), 1:max_prob_ind(i)]);
            x1_temp = [x1_straight(min_prob_ind(i):end), x1_straight(1:max_prob_ind(i))+180];
        end
        % avoid repeated values, only consider middle range
        op = 1;
        ed = length(prob_temp);
        for j = 1:length(prob_temp)
            if prob_temp(j) > 0.25
                op = j-1;
                break;
            end
        end
        for j = length(prob_temp):-1:1
            if prob_temp(j) < 0.75
                ed = j+1;
                break;
            end
        end
        %%%% check %%%%
        if op < 1
            error('Invalid psychometric starting index.')
        end
        if ed > length(x1_temp)
            error('Invalid psychometric ending index.')
        end
        C = unique(prob_temp(op:ed));
        if length(C)<ed-op+1 % there are repeated values
            error('Repeated values in psychometric curve.')
        end
        %%%%%%%%%%%%%%%
        POE1(i) = interp1(prob_temp(op:ed), x1_temp(op:ed), 0.5, 'spline', 'extrap');
        thresh0 = interp1(prob_temp(op:ed), x1_temp(op:ed), [0.25, 0.75], 'spline', 'extrap');
        thresh(i) = (thresh0(2) - thresh0(1))/2;
    end
    bias = circ90(POE1-stimulus);
end

end
