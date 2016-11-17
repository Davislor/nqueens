#define NDEBUG 1
#define _XOPEN_SOURCE 600

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>

#if defined(USE_C99_EXTENSIONS) || ( __GNUC__ >= 4 ) || ( _STDC_VERSION__ >= 199901L )
#  define USE_C99_EXTENSIONS 1
#else
#  warning "If your compiler supports C99, define USE_C99_EXTENSIONS.  Thank you."
#endif

#ifdef USE_C99_EXTENSIONS

#include <stdint.h>
#include <math.h>

#else

#include <cmath>

typedef unsigned long uint_fast32_t;
typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned long uint32_t;

using std::frexp;

inline int ilogb ( double x )
{
  int result;

  frexp( x, & result );

  return --result;
}

#endif

using std::clock;
using std::cout;
using std::endl;
using std::ostream;
using std::setw;
using std::transform;
using std::vector;

typedef uint_fast32_t uint_fast;

/* This template contains a generic implementation of a N-queens algorithm.
 * The compiler should be able to optimize for speed (at the expense of code
 * size, which is far less important) based on the fact that N is a constant.
 * For efficiency, we can plug in hand-optimized implementations for any
 * value of N.
 *
 * The algorithm assumes that an unsigned long has at least N bits.  Given
 * that it has at least 32 bits, and any computer capable of solving the
 * problem for N > 32 in any feasible amount of time will have at least a 64-
 * bit architecture, this seems safe.
 *
 * The largest problem with this program was its choice of data types; the
 * current solution is to make bitfield a template parameter.  It should be
 * an unsigned integral type, although you might see whether unsigned long
 * is faster on your architecture.
 */

template < typename T >
  class bitwise_or
  {
    // A binary functor, intended for bitfields.
    public:

    inline T operator() ( const T& a, const T& b )
    {
      return (a | b);
    }
  };

template < unsigned int N, typename bitfield >
  class board
  {
/* This object represents a N * N chessboard.  Internally, it maintains a
 * vector of bitfields, each of size N.  A modern compiler should be able to
 * vectorize operations on these data.  On most implementations, this should
 * fit each bitfield into the smallest appropriate integer type, permitting
 * a SIMD instruction such as POR.  This is probably not the absolutely most
 * compact representation (e.g., the 5x5 board will probably need 5 bytes
 * instead of 4), but it should generate appropriately-aligned code for
 * speed.
 *
 * A 1 in the position representing any square means that a queen definitely
 * threatens that square.  A 0 means that we do not know of any queen which
 * threatens that square.
 *
 * The bit positions are stored in row-major order, with the square on the
 * upper left represented by squares[0][0].
 */
    public:

    class internal
    {
      public:

      inline internal(void) : bits(0U) {};
      inline internal( bitfield x ) : bits(x) {};
      inline operator uint_fast(void) const { return (uint_fast)bits; }
      inline const internal operator| ( const internal& x ) const
        { return internal ( bitfield (bits | x.bits) ); }

      private:

      bitfield bits;

      friend class board;
    };

    typedef vector< vector< internal > > attack_matrix_t;

// We should not need to override the implicit member functions.

    void clear(void);
// Wipes all queens off the board;
    void place_queen( uint_fast row, uint_fast column_bit );
// Place a queen, updating the board as to which squares it attacks.
    uint_fast first_safe_space( uint_fast row ) const;
// Returns the first column in the given row not under attack, or
// a number out of range if no such column exists.
    uint_fast next_safe_space ( uint_fast row,
                                uint_fast column_bit
                              ) const;
/* Returns either the first column to the right of the specified column
 * that is not under attack, or a number out of range if no such column
 * exists.
 */
    static void initialize_attack_matrix(void);
    static void destroy_attack_matrix(void);

    private:

    static attack_matrix_t attack_matrix;

    internal squares[N];

// Stores N * N bits representing the positions under attack on a N * N
// chessboard.
  };

template < unsigned int N, typename bitfield >
   typename board<N, bitfield>::attack_matrix_t board< N, bitfield >::attack_matrix;

/* Stores bitfields representing positions attacked by a queen in a given
 * square.  The place_queen() member function fills the individual matrix
 * entries in as needed.
 *
 * A queen in row i, column j corresponds to the vector in
 * attack_matrix[ i*N + j ].  There are no entries for the last row, N-1,
 * because they should never be needed.  Only the bitfields representing rows
 * i+1 through N-1 are stored.
 */

template < unsigned int N, typename bitfield >
  void board< N, bitfield >::clear( void )
// Wipes all queens off the board.
{
  for ( size_t i = 0; i < N; ++i )
    squares[i].bits = 0;
}

template < unsigned int N, typename bitfield >
  void board< N, bitfield >::place_queen( uint_fast row,
                                          uint_fast column_bit
                                        )
// Place a queen on the board, updating all rows below row to show which
// squares are under attack.
  {
/* One annoyance of representing column by a bitfield is that, to convert it
 * back to an index, we use floating-point math.  This is still a constant-
 * time operation on modern computers.  Furthermore, finding the next column
 * to search is also a constant-time operation.  Therefore, this method
 * beats linear search for large N.
 *
 * We aren't using large N, so whether this is really faster is doubtful.
 */
    const uint_fast column =
static_cast<uint_fast> ( ilogb( static_cast<double> (column_bit) ) );
    const ptrdiff_t i = N * row + column;
    static const bitwise_or<internal> or_functor = bitwise_or<internal>();

// If the test run passes all these sanity checks, it should be safe to
// recompile with NDEBUG, and run without them.
    assert( row >= 0 );
    assert( row < N - 1 );
    assert( column >= 0 );
    assert( column < N );

// If the requested square is already in the database, don't enter it again.
    if ( attack_matrix[i].empty() )
    {
      vector< internal > attacks;

      attacks.reserve( static_cast<size_t>(N - 1 - row) );

// Find the highest square below us that is in the matrix, then fill upwards
// from there.
      ptrdiff_t j;
      uint_fast shift = 0;

      for ( j = i + N; j < N * N - N; j += N )
        if ( ! attack_matrix[j].empty() )
          {
            attacks = attack_matrix[j];
            shift = static_cast<uint_fast> (attacks.size());

            break;
          }

// If we fell through the loop, we are on the final row.  If we didn't, we
// are on the last filled square.  In either case, we want to move one square
// up.

/* Starting with the second-to-last row, fill in the attack matrix for the
 * square in that row directly beneath the one we want.
 */
     for ( j -= N; j >= i; j -= N )
     {
       ++shift;

       const uint_fast bottom_row ( column_bit |
                                    column_bit << shift |
                                    column_bit >> shift
                                  );

       attacks.push_back( internal(static_cast<bitfield>(bottom_row) ) );
       attack_matrix[j] = attacks;
     }

// Now that we have the attack matrix, AND it with the board.

     transform ( (&squares[0]) + 1 + row,
                 (&squares[0]) + N,
                 &attacks[0],
                 (&squares[0]) + 1 + row,
                 or_functor
               );

      assert ( attacks.size() == N - 1 - row; );
    } // end if ( attack_matrix[i].empty() )
    else
    {
     transform ( (&squares[0]) + 1 + row,
                 (&squares[0]) + N,
                 &(attack_matrix[i][0]),
                 (&squares[0]) + 1 + row,
                 or_functor
               );
    }
  }

template < unsigned int N, typename bitfield >
  uint_fast board< N, bitfield >::first_safe_space( uint_fast row ) const
  {
    const uint_fast to_return = squares[row];

    return ( ~to_return & (to_return + 1) );
  }

template < unsigned int N, typename bitfield >
  uint_fast board< N, bitfield >::next_safe_space( uint_fast row,
                                                   uint_fast column_bit
                                                 ) const
  {
    uint_fast to_return = squares[row];

// Mask the columns up to column_bit.
    to_return |= column_bit;
    to_return |= column_bit - 1;

    return ( ~to_return & (to_return + 1) );
  }


template < unsigned int N, typename bitfield >
  void board< N, bitfield >::initialize_attack_matrix (void)
  {
    attack_matrix.resize( (N - 1) * N );
  }

template < unsigned int N, typename bitfield >
  void board< N, bitfield >::destroy_attack_matrix (void)
  {
    attack_matrix.clear();
  }
 

template < unsigned int N, typename bitfield >
  long solve_queens(void)
  {
    long solutions = 0;
    uint_fast j[N-1];
    board< N, bitfield > backtrack[N-1];

    board< N, bitfield >::initialize_attack_matrix();

/* Find the solutions where the queen on the first row is on the left half
 * of the board.  For each such solution, there is a one-to-one corresp-
 * ondence with the solutions where that same queen is on the right half
 * of the board, and the solutions are reflections of each other.
 *
 * If N is odd, any solution with the last queen in the center column will
 * have at least three counterparts: its horizontal reflection, its vertical
 * reflection, and its horizontal and vertical reflection.  Any solution with
 * the first queen in that column would be one of the last two.  Therefore,
 * we can eliminate the loop to test the center column.
 */
    for ( j[0] = 1; j[0] < 1 << (N >> 1); j[0] <<= 1 )
    {
      backtrack[0].clear();
      backtrack[0].place_queen( 0, j[0] );

      uint_fast current_row = 1;
      j[1] = backtrack[0].first_safe_space(1);

      while (true)
      {
        if ( 1 << N <= j[current_row] )
        {
// There are no (more) possible matches on this row.
          if ( 1 == current_row )
            break;

          --current_row;

          j[current_row] =
backtrack[current_row-1].next_safe_space( current_row, j[current_row] );
        }
        else if ( N - 2 == current_row )
        {
// Descend to the final level.
//
// Any open spaces we find here (there should be one at most) is a
// solution.  Furthermore, it has a symmetric solution.
          backtrack[N-2] = backtrack[N-3];
          backtrack[N-2].place_queen( N-2, j[N-2] );

          const uint_fast k = backtrack[N-2].first_safe_space(N-1);

          if ( k < 1 << N )
        {
          if ( N % 2 == 1 && k == 1 << (N >> 1) )
// Our last queen is in the center column.
            solutions += 4;
          else
            solutions += 2;
        }

          j[N-2] = backtrack[N-3].next_safe_space( N-2, j[N-2] );
        }
        else
        {
// Descend.
          backtrack[current_row] = backtrack[current_row - 1];
          backtrack[current_row].place_queen(current_row, j[current_row] );

          ++current_row;
          j[current_row] =
backtrack[current_row-1].first_safe_space(current_row);
        } // end if
      } // end while
    } // end for

    board< N, bitfield >::destroy_attack_matrix();

    return solutions;
  }

// For the degenerate case N <= 2, we need a trivial handler:
template < >
  long solve_queens< 1, uint8_t > (void)
  {
    return 1;
  }

template < >
  long solve_queens< 2, uint8_t > (void)
  {
    return 0;
  }

long queens( unsigned n )
{
  switch (n)
  {
    case  1: return solve_queens<  1, uint8_t  > (); break;
    case  2: return solve_queens<  2, uint8_t  > (); break;
    case  3: return solve_queens<  3, uint8_t  > (); break;
    case  4: return solve_queens<  4, uint8_t  > (); break;
    case  5: return solve_queens<  5, uint8_t  > (); break;
    case  6: return solve_queens<  6, uint8_t  > (); break;
    case  7: return solve_queens<  7, uint8_t  > (); break;
    case  8: return solve_queens<  8, uint8_t  > (); break;
    case  9: return solve_queens<  9, uint16_t > (); break;
    case 10: return solve_queens< 10, uint16_t > (); break;
    case 11: return solve_queens< 11, uint16_t > (); break;
    case 12: return solve_queens< 12, uint16_t > (); break;
    case 13: return solve_queens< 13, uint16_t > (); break;
    case 14: return solve_queens< 14, uint16_t > (); break;
    case 15: return solve_queens< 15, uint16_t > (); break;
    case 16: return solve_queens< 16, uint16_t > (); break;
    case 17: return solve_queens< 17, uint32_t > (); break;
    case 18: return solve_queens< 18, uint32_t > (); break;

    default: return -1;
// If you're getting -1, expand the switch block.
  }
}

int main(void)
{
  static const unsigned start = 1;
//  static const unsigned finish = 14;
  static const unsigned finish = 17;

  clock_t begin, end;
  clock_t total = 0;

  unsigned i;

  for ( i = start; i <= finish; ++i )
  {
    long answer;

    begin = clock();
    answer = queens(i);
    end = clock();

    total += end - begin;

    cout << setw(8);
    cout << static_cast<double>( 1000 * (end - begin) ) / CLOCKS_PER_SEC;
    cout << " ms # of solutions for ";
    cout << setw(2) << i << " x ";
    cout << setw(2) << i << " = ";
    cout << answer << endl;
  }

  cout << "Total time: ";
  cout << static_cast<double>( 1000 * total ) / CLOCKS_PER_SEC;
  cout << " ms." << endl;

  return EXIT_SUCCESS;
}
