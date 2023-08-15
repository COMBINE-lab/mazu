pub type Result<T> = std::result::Result<T, Error>;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("failed to Elias-Fano encode not-monotone sequence")]
    EFNotMonotone,

    #[error("cannot Elias-Fano encode empty sequence")]
    EFEmpty,

    #[error("IO error")]
    IOError(#[from] std::io::Error),

    #[error("Serde JSON error")]
    SerdeJSON(#[from] serde_json::Error),

    #[error("Bincode error")]
    Bincode(#[from] bincode::Error),

    #[error("failed to load index")]
    IndexLoad,

    #[error("invalid data: {0}")]
    InvalidData(String),

    #[error("failed with error: {0}")]
    Other(String),

    #[error("failed to parse cuttlefish tiling token")]
    CfSeqTokenParseError,
}

impl Error {
    pub fn other<T: ToString>(msg: T) -> Self {
        Self::Other(msg.to_string())
    }
}
