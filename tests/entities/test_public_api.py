import protein_graph.entities as entities


def test_all_declared_names_are_importable() -> None:
    for name in entities.__all__:
        assert hasattr(entities, name), (
            f"{name} declared in __all__ but not importable"
        )
