****************************
Design of Transformed Module
****************************

Some quick notes about the design of the ``Transformed`` module.

- ``property_types::Transformed<T>` is templated on the property type ``T``
   
   - This allows modules to specify they need a transformed overlap matrix by``Transformed<Overlap>``