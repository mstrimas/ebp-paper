
validate <- function(model, data) {
  # predict on test data set
  if (inherits(model$model, "maxnet")) {
    print("predict maxnet")
    # predict on standardised fixed test dataset
    pred <- predict(model$model, newdata = data, 
                    type = "logistic", clamp = FALSE)
    pred <- data.frame(id = nrow(data), 
                       obs = data$species_observed,
                       pred = pred)
  } else {
    pred_rf <- predict(model$model, data = data, type = "response")
    pred_cal <- predict(model$calibration, 
                        newdata = data.frame(pred = pred_rf$predictions[,2]), 
                        type = "response")
    pred_cal <- ifelse(ifelse(pred_cal<0, 0, pred_cal)>1, 1, pred_cal)
    pred <- data.frame(id = nrow(data),
                       obs = data$species_observed,
                       pred = pred_cal)
  }
  pred <- drop_na(pred)
  
  # validation metrics
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
         specificity = round(pa_metrics$specificity, 5)
  )
}

