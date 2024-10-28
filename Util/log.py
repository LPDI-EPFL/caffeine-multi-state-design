def log(logFile, string):
  logFile.write(string+"\n")
  print string
  logFile.flush()
