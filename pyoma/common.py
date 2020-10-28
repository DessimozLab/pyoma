import logging

package_logger = logging.getLogger("pyoma")
package_logger.addHandler(logging.NullHandler())
