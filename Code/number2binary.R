#from https://stackoverflow.com/questions/12088080/how-to-convert-integer-number-into-binary-vector

number2binary = function(number, noBits) {
       binary_vector = rev(as.numeric(intToBits(number)))
       if(missing(noBits)) {
          return(binary_vector)
       } else {
          binary_vector[-(1:(length(binary_vector) - noBits))]
       }
    }
