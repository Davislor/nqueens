#define NDEBUG 1

#include <algorithm>
#include <bitset>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

#include <math.h>

using std::bitset;
using std::clock;
using std::cout;
using std::endl;
using std::ostream;
using std::setw;
using std::transform;
using std::vector;

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
 */

typedef unsigned long column_t;
// A column_t is a bitfield with at most a single bit set, in the bit
// position corresponding to the specified column.

// This functor performs a bitwise AND on its two predicates.
template <unsigned int N>
  struct bitwise_or
  {
    public:

    inline bitset<N> operator() ( const bitset<N> a, const bitset<N> b );
  };

template <unsigned int N>
  inline bitset<N> bitwise_or<N>::operator() ( const bitset<N> a,
                                               const bitset<N> b
                                             )
  {
    return bitset<N> ( a | b );
  }

template <unsigned int N>
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
// We should not need to override the implicit member functions.

    void clear (void);
// Wipes all queens off the board;
    void place_queen( unsigned row, column_t column_bit );
// Place a queen, updating the board as to which squares it attacks.
    column_t first_safe_space( unsigned row ) const;
// Returns the first column in the given row not under attack, or
// (1 << N) if no such column exists.
    column_t next_safe_space ( unsigned row,
                               column_t column_bit
                             ) const;
/* Returns either the first column to the right of the specified column
 * that is not under attack, or (1 << N) if no such column exists.
 */
    static void initialize_attack_matrix(void);
    static void destroy_attack_matrix(void);

    private:
    static vector< vector< bitset<N> > > attack_matrix;

    bitset<N> squares[N];
// Stores N * N bits representing the positions under attack on a N * N
// chessboard.
  };

template <unsigned int N>
  vector< vector< bitset<N> > > board<N>::attack_matrix;
/* Stores bitfields representing positions attacked by a queen in a given
 * square.  The place_queen() member function fills the individual matrix
 * entries in as needed.
 *
 * A queen in row i, column j corresponds to the vector in
 * attack_matrix[ i*N + j ].  There are no entries for the last row, N-1,
 * because they should never be needed.  Only the bitfields representing rows
 * i+1 through N-1 are stored.
 */

template <unsigned int N>
  void board<N>::clear( void )
// Wipes all queens off the board.
{
  for ( size_t i = 0; i < N; ++i )
    squares[i].reset();
}

template <unsigned int N>
  void board<N>::place_queen( unsigned row, column_t column_bit )
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
    const unsigned column =
static_cast<unsigned> ( ilogb( static_cast<double> (column_bit) ) );
    const ptrdiff_t i = N * row + column;
    static const bitwise_or<N> or_functor = bitwise_or<N>();

// If the test run passes all these sanity checks, it should be safe to
// recompile with NDEBUG, and run without them.
    assert( row >= 0 );
    assert( row < N - 1 );
    assert( column >= 0 );
    assert( column < N );
    assert( squares.size() == N );
#ifndef NDEBUG
    {
      const bitset<sizeof(unsigned long) * CHAR_BIT> temp(column_bit);
      assert ( 1 == temp.count() );
    }
#endif

    const ptrdiff_t squares_to_fill = N - 1 - row;

// If the requested square is already in the database, don't enter it again.
    if ( attack_matrix[i].empty() )
    {
// Find the highest square below us that is in the matrix, then fill upwards
// from there.
      ptrdiff_t j;
      unsigned long shift = 0;
      vector< bitset<N> > attacks;

      attacks.reserve( static_cast<size_t> (squares_to_fill) );

      for ( j = i + N; j < N * N - N; j += N )
        if ( ! attack_matrix[j].empty() )
          {
            attacks = attack_matrix[j];
            shift = attacks.size();

            break;
          }

// If we fell through the loop, we are on the final row.  If we didn't, we
// are on the last filled square.  In either case, we want to move one square
// up.

      j -= N;

      {
/* Starting with the second-to-last row, fill in the attack matrix for the
 * square in that row directly beneath the one we want.
 */
        for ( ; j >= i; j -= N )
        {
          ++shift;

          const bitset<N> bottom_row = column_bit |
                                        column_bit << shift |
                                        column_bit >> shift;

          attacks.push_back(bottom_row);
          attack_matrix[j] = attacks;
        } // end for
      } // end if

      assert ( attacks.size() == (size_t) squares_to_fill );
// Now that we have the attack matrix, AND it with the board.

      transform ( &squares[0] + 1 + row,
                  &squares[0] + N,
                  &attacks[0],
                  &squares[0] + 1 + row,
                  or_functor
                );
    } // end if ( attack_matrix[i].empty() )
    else
    {
      transform ( &squares[0] + 1 + row,
                  &squares[0] + N,
                  &(attack_matrix[i][0]),
                  &squares[0] + 1 + row,
                  or_functor
                );
    }
  }

template <unsigned int N>
  column_t board<N>::first_safe_space( unsigned row ) const
  {
    const unsigned long to_return = squares[row].to_ulong();

    return static_cast<column_t> ( ~to_return & (to_return + 1) );
  }

template <unsigned int N>
  column_t board<N>::next_safe_space( unsigned row, column_t column_bit ) const
  {
    unsigned long to_return = squares[row].to_ulong();

#ifndef NDEBUG
    {
      bitset<sizeof(column_t) * CHAR_BIT> temp(column_bit);
      assert ( 1 == temp.count() );
    }
#endif

// Mask the columns up to column_bit.
    to_return |= column_bit;
    to_return |= column_bit - 1;

    return static_cast<column_t> ( ~to_return & (to_return + 1) );
  }

template <unsigned int N>
  void board<N>::initialize_attack_matrix (void)
  {
    attack_matrix.resize( (N - 1) * N );
  }

template <unsigned int N>
  void board<N>::destroy_attack_matrix (void)
  {
    attack_matrix.clear();
  }

#ifndef NDEBUG
template <unsigned int N>
  void print_board ( const board<N> &x )
{
  for ( size_t i = 0; i < N; ++i )
  {
    column_t t = x.first_safe_space(i);

    for ( column_t j = 1; j < 1 << N; j <<=1 )
    {
      if ( t == j )
      {
        cout << 'O';
        t = x.next_safe_space(i, j);
      }
      else
        cout << 'X';
    }
    cout << endl;
  }
  cout << endl;

  return ;
}
#endif

template <unsigned int N>
  long solve_queens(void)
  {
    long solutions = 0;
    column_t j[N-1];
    board<N> backtrack[N-1];

    board<N>::initialize_attack_matrix();

/* Find the solutions where the queen on the first row is on the left half
 * of the board.  For each such solution, there is a one-to-one corresp-
 * ondence with the solutions where that same queen is on the right half
 * of the board, and the solutions are reflections of each other.
 */
    for ( j[0] = 1; j[0] < 1 << (N/2); j[0] <<= 1 )
    {
      backtrack[0].clear();
      backtrack[0].place_queen( 0, j[0] );

      size_t current_row = 1;
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

          if ( backtrack[N-2].first_safe_space(N-1) < 1 << N )
            solutions += 2;

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

// If N is odd, we also need to check the middle column.  The two halves of
// the board on either side are symmetric, so we need only check one.
    if ( 1 == N % 2 )
    {
#ifndef NDEBUG
      j[0] = 1 << (N/2);
#endif

      backtrack[0].clear();
      backtrack[0].place_queen( 0, 1 << (N/2) );

// The second queen cannot be in the position diagonally below the first.
      for ( j[1] = 1; j[1] < 1 << (N/2 - 1); j[1] <<= 1 )
      {
        backtrack[1] = backtrack[0];
        backtrack[1].place_queen( 1, j[1] );

        size_t current_row = 2;
        j[2] = backtrack[1].first_safe_space(2);

        while (true)
        {
          if ( 1 << N <= j[current_row] )
          {
// There are no (more) possible matches on this row.
            if ( 2 == current_row )
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

            if ( backtrack[N-2].first_safe_space(N-1) < 1 << N )
              solutions += 2;

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
    } // end if

    board<N>::destroy_attack_matrix();

    return solutions;
  }

// For the degenerate cases N < 4, we need trivial handlers:
template <>
  long solve_queens<1> (void)
  {
    return 1;
  }

long queens( unsigned n )
{
  switch (n)
  {
    case  1: return solve_queens<1> (); break;
    case  2: return solve_queens<2> (); break;
    case  3: return solve_queens<3> (); break;
    case  4: return solve_queens<4> (); break;
    case  5: return solve_queens<5> (); break;
    case  6: return solve_queens<6> (); break;
    case  7: return solve_queens<7> (); break;
    case  8: return solve_queens<8> (); break;
    case  9: return solve_queens<9> (); break;
    case 10: return solve_queens<10>(); break;
    case 11: return solve_queens<11>(); break;
    case 12: return solve_queens<12>(); break;
    case 13: return solve_queens<13>(); break;
    case 14: return solve_queens<14>(); break;
    case 15: return solve_queens<15>(); break;
    case 16: return solve_queens<16>(); break;

    default: return -1;
// If you're getting -1, expand the switch block.
  }
}

int main(void)
{
  static const unsigned start = 1;
  static const unsigned finish = 16;
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
  cout << "ms." << endl;

  return EXIT_SUCCESS;
}
