# print and summary output

    Code
      print(summary(fit), digits = 4)
    Output
      
      Call:
      serial_fgls(formula = y ~ x, data = sim)
      
      Method:               FGLS (serial)
      Covariance:           robust
      Observations:         200 (effective 199)
      Inference:            normal
      AR order k*:          1 [BIC, kmax = 12]
      AR coefficients:      0.3825
      Non-exogeneity share: 0.035
      
      Coefficients:
                  Estimate Std. Error z value Pr(>|z|)    
      (Intercept)  0.95364    0.10739    8.88   <2e-16 ***
      x            1.08684    0.05719   19.00   <2e-16 ***
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      Notes:
      [1] Covariance: PG26 eq. (29), robust to non-exogenous (pre-determined)
          regressors.
      [2] Collinear lagged columns dropped from the Durbin regression:
          L1.(Intercept).
      

---

    Code
      print(fit, digits = 4)
    Output
      
      Call:
      serial_fgls(formula = y ~ x, data = sim)
      
      FGLS (serial); AR order k* = 1 (BIC), rho = 0.3825; covariance: robust
      
      Coefficients:
      (Intercept)            x  
           0.9536       1.0868  
      

