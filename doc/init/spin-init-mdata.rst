dpgen spin_init machine parameters
==================================

.. dargs::
   :module: dpgen.data.arginfo
   :func: spin_init_mdata_arginfo

The generated schema above describes the ``fp``/``spin`` VASP entries.
A Stage-5-only run does not require ``fp``; it requires a top-level
``convert-data`` entry in MACHINE.
It accepts the same ``command`` / ``machine`` / ``resources`` fields as an
``fp`` entry, either directly as an object or inside a one-element list. The
converter must produce ``out/data.extxyz`` from each scale's ``data/`` input.
See :doc:`spin-init-usage` for a complete example.
