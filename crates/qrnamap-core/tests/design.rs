//! Design formula grammar and contrast tests.

use qrnamap_core::{ColType, Contrast, Design, DesignError};

fn cols() -> Vec<(String, ColType)> {
    vec![
        ("condition".to_string(), ColType::Factor),
        ("batch".to_string(), ColType::Factor),
        ("dose".to_string(), ColType::Numeric),
    ]
}

fn bad_syntax(f: &str) {
    match Design::parse(f, &cols()) {
        Err(DesignError::BadSyntax(_)) => {}
        other => panic!("{f:?} gave {other:?}"),
    }
}

#[test]
fn accepts_condition() {
    let d = Design::parse("~ condition", &cols()).unwrap();
    assert_eq!(d.to_r(), "~ condition");
    assert_eq!(d.terms(), ["condition"]);
}

#[test]
fn accepts_batch_plus_condition() {
    let d = Design::parse("~ batch + condition", &cols()).unwrap();
    assert_eq!(d.to_r(), "~ batch + condition");
    assert_eq!(d.terms(), ["batch", "condition"]);
}

#[test]
fn whitespace_is_ignored() {
    let d = Design::parse("  ~batch+   condition ", &cols()).unwrap();
    assert_eq!(d.to_r(), "~ batch + condition");
    let d = Design::parse("~\tdose", &cols()).unwrap();
    assert_eq!(d.to_r(), "~ dose");
}

#[test]
fn interaction_and_star_terms_accepted() {
    let d = Design::parse("~ batch + condition:dose", &cols()).unwrap();
    assert_eq!(d.to_r(), "~ batch + condition:dose");
    let d = Design::parse("~ condition * batch", &cols()).unwrap();
    assert_eq!(d.to_r(), "~ condition*batch");
}

#[test]
fn rejects_unknown_column() {
    assert_eq!(
        Design::parse("~ treatment", &cols()),
        Err(DesignError::UnknownColumn("treatment".into()))
    );
    assert_eq!(
        Design::parse("~ batch + nope", &cols()),
        Err(DesignError::UnknownColumn("nope".into()))
    );
    // Column names are case sensitive.
    assert_eq!(
        Design::parse("~ Condition", &cols()),
        Err(DesignError::UnknownColumn("Condition".into()))
    );
}

#[test]
fn rejects_function_call() {
    bad_syntax("~ log(x)");
    bad_syntax("~ log(dose)");
    bad_syntax("~ condition + factor(batch)");
}

#[test]
fn rejects_semicolon_and_backtick() {
    bad_syntax("~ condition; system('id')");
    bad_syntax("~ `condition`");
    bad_syntax("~ condition`");
}

#[test]
fn rejects_other_bad_input() {
    for f in [
        "",
        "~",
        "condition",
        "condition ~ batch",
        "~ ~ condition",
        "~ condition +",
        "~ + condition",
        "~ condition + + batch",
        "~ condition batch",
        "~ condition - batch",
        "~ condition | batch",
        "~ $condition",
        "~ \"condition\"",
        "~ 'condition'",
        ".",
        "~ .",
        "~ 1",
        "~ condition + 0",
        "~ condition\nbatch",
        "~ cond\\ition",
        "~ (condition)",
        "~ _x",
        "~ condition # c",
    ] {
        bad_syntax(f);
    }
}

#[test]
fn syntax_error_does_not_echo_long_user_text() {
    let long = "~ condition + ".to_string() + &"$".repeat(500);
    match Design::parse(&long, &cols()) {
        Err(DesignError::BadSyntax(m)) => {
            assert!(m.len() < 100, "message too long: {m}");
            assert!(!m.contains("$$$"));
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn to_r_quotes_nothing_but_names_and_operators() {
    let c = vec![
        ("condition".to_string(), ColType::Factor),
        ("batch".to_string(), ColType::Factor),
        ("dose".to_string(), ColType::Numeric),
        ("a".to_string(), ColType::Factor),
        ("b".to_string(), ColType::Factor),
    ];
    for f in [
        "~ condition",
        "~ batch + condition",
        "~ batch+condition:dose",
        "~ a*b",
    ] {
        let r = Design::parse(f, &c).unwrap().to_r();
        assert!(r.starts_with('~'));
        assert!(
            r.chars()
                .all(|ch| ch.is_ascii_alphanumeric() || " _.~+:*".contains(ch)),
            "{r}"
        );
        assert!(!r.contains('"') && !r.contains('\'') && !r.contains('`'));
    }
    // Dots and underscores in names pass through unchanged.
    let c = vec![("b_1.x".to_string(), ColType::Factor)];
    assert_eq!(Design::parse("~ b_1.x", &c).unwrap().to_r(), "~ b_1.x");
}

#[test]
fn contrast_parse_valid_and_invalid() {
    let c = Contrast::parse("condition,treat,ctrl", &cols()).unwrap();
    assert_eq!(
        c,
        Contrast {
            factor: "condition".into(),
            level: "treat".into(),
            reference: "ctrl".into()
        }
    );
    let c = Contrast::parse(" condition , t-1.a_b , c0 ", &cols()).unwrap();
    assert_eq!(c.level, "t-1.a_b");
    assert_eq!(c.reference, "c0");

    assert_eq!(
        Contrast::parse("nope,a,b", &cols()),
        Err(DesignError::UnknownColumn("nope".into()))
    );
    // Numeric column cannot be the factor.
    assert!(matches!(
        Contrast::parse("dose,a,b", &cols()),
        Err(DesignError::BadSyntax(_))
    ));
    for bad in [
        "condition,tr\"eat,ctrl",
        "condition,treat,ct$rl",
        "condition,tr eat,ctrl",
        "condition,treat,ct`rl",
        "condition,tr\\eat,ctrl",
        "condition,'x',ctrl",
        "condition,treat",
        "condition,treat,ctrl,extra",
        "condition,,ctrl",
        "condition,treat,",
        "condition,same,same",
        "",
        "cond`ition,a,b",
    ] {
        assert!(
            matches!(
                Contrast::parse(bad, &cols()),
                Err(DesignError::BadSyntax(_))
            ),
            "{bad:?} -> {:?}",
            Contrast::parse(bad, &cols())
        );
    }
}

#[test]
fn errors_display() {
    assert!(DesignError::UnknownColumn("x".into())
        .to_string()
        .contains('x'));
    let _: &dyn std::error::Error = &DesignError::BadSyntax("m".into());
}
