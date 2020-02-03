class Helpers {

  def containsValueRecursive(def top, query) {
      for(it in top) {
          if (it.value instanceof Map) {
              if(containsValueRecursive(it.value, query))
                return true
          } else if (it.value == query) {
              return true
          }
      }
      return false
  }

}