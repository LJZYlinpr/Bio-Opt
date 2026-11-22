"""Stable filenames and random seeds for physical-depth targets."""


def target_layers(pair):
    depths = pair["targets"]
    if "target_indices" not in pair:
        return [(str(z), z, 42 + int(z // 20)) for z in depths]
    return [(f"layer_{i:04d}", z, 42 + i)
            for i, z in zip(pair["target_indices"], depths, strict=True)]
