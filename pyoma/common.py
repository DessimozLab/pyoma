import logging

DEBUG_LOG_OBJECT_INDENT_LEN = 4


package_logger = logging.getLogger('pyoma')
package_logger.addHandler(logging.NullHandler())

debug_log_object_indent = ' ' * DEBUG_LOG_OBJECT_INDENT_LEN

def debug_log_object_dump(obj):
    return re.sub(r'(^|\n)', r'\1' + debug_log_object_indent,
                      pprint.pformat(obj, width=(80 - DEBUG_LOG_OBJECT_INDENT_LEN)))



