  
validate <- function(model, data, bbs_combine = FALSE) {

  pred_vec <- predict_bp_model(model, data)
  pred <- data.frame(id = nrow(data), 
                     obs = data$species_observed,
                     pred = pred_vec)

  pred <- drop_na(pred)

  if(bbs_combine){
    
  }
  
  # validation metrics
  print("calculate ppms")
  # mean squared error (mse)
  mse <- mean((pred$obs - pred$pred)^2, na.rm = TRUE)
  # pick threshold to maximize kappa
  thresh <- optimal.thresholds(pred, opt.methods = "MaxKappa")[1, 2]
  # calculate accuracy metrics: auc, kappa, sensitivity, specificity, brier
  pa_metrics <- presence.absence.accuracy(pred, threshold = thresh, 
                                          na.rm = TRUE, st.dev = FALSE)
  
  # summarise the performance of this model
  tibble(mse = round(mse, 5),
         auc = round(pa_metrics$AUC, 5),
         # threshold required
         threshold = round(thresh, 5),
         kappa = round(pa_metrics$Kappa, 5), 
         sensitivity = round(pa_metrics$sensitivity, 5),
         specificity = round(pa_metrics$specificity, 5),
         n_checklists = model$n_checklists, 
         n_pos = model$n_sightings
  )
}
