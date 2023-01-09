pub enum BaseMapError {
    KeyError,
    BaseError,
    IndexError,
    IntegerOverflowError,
    PositionError,
    BaseQualityError,
    MappingQualityError,
}

use std::{
    fmt::{Display, Formatter, Result as FmtResult},
    num::ParseIntError,
};

/// The Error definition.  Variants can carry payloads, which I am using
/// here to carry the source or causal error.  The variant does not need to
/// have the same name as the foreign error type it wraps.
#[derive(Debug)]
pub enum Error {
    // Other error variant definitions ...
    ParseIntError(ParseIntError),
}

/// Allows your error to be displayed using `{}`, and not just `{:?}`
impl Display for Error {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(
            f,
            "{}",
            match self {
                Error::ParseIntError(err) => format!("Error parsing data as an integer: {:?}", err),
                // Add display implementations for other Error variants here
            }
        )
    }
}

/// Automatically converts a `std::string::ParseError` into a
/// `<this_crate>::error::Error` whenever coercion is needed--e.g. when
/// using the unwrap/early-return operator (`?`)
impl From<ParseIntError> for Error {
    fn from(err: ParseIntError) -> Self {
        Error::ParseIntError(err)
    }
}
