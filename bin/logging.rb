module Logging
  # found this on http://stackoverflow.com/questions/917566/ruby-share-logger-instance-among-module-classes
  @out = STDERR
  @level = Logger::INFO

  def logger
    @logger ||= Logging.logger_for(self.class.name)
  end

  def configure(options)
    @out = options['logout']
    @level = options['log_level']
  end

  @loggers = {}

  class << self
    def logger_for(classname)
      @loggers[classname] ||= configure_logger_for(classname)
    end

    def configure_logger_for(classname)
      logger = Logger.new(@out)
      logger.progname = classname
      logger.level = @level
      logger
    end
  end
end