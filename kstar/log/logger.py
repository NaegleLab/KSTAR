import logging
def get_logger(name, filename):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(filename)
    log_format = logging.Formatter('%(asctime)s\t%(name)s\t%(levelname)s:\t%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    handler.setFormatter(log_format)
#     if (logger.hasHandlers()):
#         logger.handlers.clear()
    logger.addHandler(handler)
    return logger