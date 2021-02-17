The Accelerator class
======================

The :any:`Accelerator` class is top-level class in KENV. It contains all
the accelerator data, and has the high-level method :any:`compile`
that performs the accelerator.

.. autoclass:: kenv.accelerator.Accelerator
   :members: compile, add_solenoid, add_accel, add_quadrupole, add_corrector_x, add_corrector_y
