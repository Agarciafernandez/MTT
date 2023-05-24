function log_pdf_clutter=Calculate_log_clutter_pdf_nb(Nc,r_nb,p_nb,Area_tot)
%Evaluates the multi-target clutter density for negative binomial clutter

% Author: Angel F. Garcia-Fernandez

%Log Probability mass function of the cardinality (see nbinpdf.m and
%https://uk.mathworks.com/help/stats/negative-binomial-distribution.html)
%pmf_card_c_log is equivalent to log(nbinpdf(Nc,r_nb,p_nb))

pmf_card_c_log=gammaln(r_nb + Nc) - gammaln(Nc + 1) - gammaln(r_nb) ...
    + r_nb.*log(p_nb) + Nc*log1p(-p_nb);


%One option would be to calculate
%log_pdf_clutter=pmf_card_c_log+log (factorial(Nc))-Nc*log(Area_tot);
%However, we calculate the log_factorial of N_c as in the
%line below because log(factorial(Nc)) returns Inf value when Nc>170. To calculate this with
%higher accuracy, we use the relation between gamma function and
%factorial.

log_factorial=gammaln(Nc+1);
log_pdf_clutter=pmf_card_c_log+log_factorial-Nc*log(Area_tot);

end