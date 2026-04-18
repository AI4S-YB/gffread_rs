use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompatError {
    pub message: String,
    pub exit_code: i32,
    pub show_usage: bool,
}

impl CompatError {
    pub fn new(message: impl Into<String>, exit_code: i32) -> Self {
        Self {
            message: message.into(),
            exit_code,
            show_usage: false,
        }
    }

    pub fn with_usage(message: impl Into<String>, exit_code: i32) -> Self {
        Self {
            message: message.into(),
            exit_code,
            show_usage: true,
        }
    }
}

impl fmt::Display for CompatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for CompatError {}
