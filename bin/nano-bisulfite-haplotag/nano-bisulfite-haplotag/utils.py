"""
Utility functions for nano-bisulfite-haplotag
"""

import os
import logging
import argparse
from typing import Union, Optional


def str2bool(v: Union[str, bool]) -> bool:
    """Convert string to boolean value."""
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def setup_logging(log_file: str, debug: bool = False) -> None:
    """Setup logging configuration."""
    log_level = logging.DEBUG if debug else logging.INFO
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Setup file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    
    # Setup console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    
    # Setup root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)


def validate_file_exists(file_path: str, file_type: str = "file") -> None:
    """Validate that a file exists."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{file_type} not found: {file_path}")
    
    if not os.path.isfile(file_path):
        raise ValueError(f"Path is not a file: {file_path}")


def create_output_directory(output_dir: str) -> None:
    """Create output directory if it doesn't exist."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info(f"Created output directory: {output_dir}")


def get_default_snp_file() -> Optional[str]:
    """Get default SNP file path if it exists."""
    default_path = "/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/xuexiaohui/database/mm10/GATK/all_DBA_2J_SNPs_C57BL_6NJ_reference.based_on_GRCm38.sorted.txt"
    
    if os.path.exists(default_path):
        return default_path
    return None