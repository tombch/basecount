use std::error::Error;
use std::fmt::{self, Display};

#[derive(Debug)]
pub enum BaseMapError {
    KeyNotFound,
    InvalidBase,
    IntegerOverflow,
    AlignmentStartNotFound,
    AlignmentEndNotFound,
    MappingQualityNotFound,
    QualityScoreNotFound,
    RegionNameNotFound,
    ReferenceSequenceIDNotFound,
}

impl Display for BaseMapError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            BaseMapError::KeyNotFound => f.write_str("KeyNotFound"),
            BaseMapError::InvalidBase => f.write_str("InvalidBase"),
            BaseMapError::IntegerOverflow => f.write_str("IntegerOverlow"),
            BaseMapError::AlignmentStartNotFound => f.write_str("AlignmentStartNotFound"),
            BaseMapError::AlignmentEndNotFound => f.write_str("AlignmentEndNotFound"),
            BaseMapError::MappingQualityNotFound => f.write_str("MappingQualityNotFound"),
            BaseMapError::QualityScoreNotFound => f.write_str("QualityScoreNotFound"),
            BaseMapError::RegionNameNotFound => f.write_str("RegionNameNotFound"),
            BaseMapError::ReferenceSequenceIDNotFound => f.write_str("ReferenceSequenceIDNotFound"),
        }
    }
}

impl Error for BaseMapError {
    fn description(&self) -> &str {
        match *self {
            BaseMapError::KeyNotFound => "Key not found",
            BaseMapError::InvalidBase => "Invalid base",
            BaseMapError::IntegerOverflow => "Integer overflow",
            BaseMapError::AlignmentStartNotFound => "Alignment start not found",
            BaseMapError::AlignmentEndNotFound => "Alignment end not found",
            BaseMapError::MappingQualityNotFound => "Mapping quality not found",
            BaseMapError::QualityScoreNotFound => "Quality score not found",
            BaseMapError::RegionNameNotFound => "Region name not found",
            BaseMapError::ReferenceSequenceIDNotFound => "Reference sequence ID not found",
        }
    }
}
